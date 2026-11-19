"""Exact measured-genotype graph calculations, reused from ANA-068."""
from pathlib import Path
import ctypes,hashlib,json,subprocess,sys,tempfile,functools,os
import numpy as np
from scipy.stats import rankdata
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components
HERE=Path(__file__).resolve().parent
STEPS=128
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def dump(p,o):Path(p).write_text(json.dumps(o,indent=2,allow_nan=False)+'\n')

def kernel():
    source=HERE/'landscape_graph_kernel.cpp';key=sha(source)[:16]
    lib=Path(tempfile.gettempdir())/('adaptenv-graph-'+key+('.dylib' if sys.platform=='darwin' else '.so'))
    if not lib.exists():
        compiler=os.environ.get('ADAPTENV_GRAPH_COMPILER')
        cmd=([compiler,'c++'] if compiler else ['clang++' if sys.platform=='darwin' else 'g++'])+['-O3','-std=c++17','-shared','-fPIC',str(source),'-o',str(lib)]
        if sys.platform=='darwin':
            sdk=subprocess.check_output(['xcrun','--show-sdk-path'],text=True).strip()
            cmd+=['-isysroot',sdk,'-isystem',str(Path(sdk)/'usr/include/c++/v1')]
        result=subprocess.run(cmd,capture_output=True,text=True)
        if result.returncode:raise RuntimeError('Graph-kernel compilation failed: '+result.stderr)
    dll=ctypes.CDLL(str(lib));D=np.ctypeslib.ndpointer(dtype=np.float64,flags='C_CONTIGUOUS');I=np.ctypeslib.ndpointer(dtype=np.int32,flags='C_CONTIGUOUS');B=np.ctypeslib.ndpointer(dtype=np.uint8,flags='C_CONTIGUOUS')
    dll.analyse.argtypes=[ctypes.c_int]*3+[D,B]+[ctypes.c_double]*3+[D,I,I]
    dll.walks.argtypes=[ctypes.c_int]*3+[D,B,I]+[ctypes.c_int]*2+[ctypes.c_double,ctypes.c_int,ctypes.c_uint64,I]
    return dll
KERNEL=None

def graph(f,valid,eps,top,floor,q=20,dim=4):
    global KERNEL
    if KERNEL is None:KERNEL=kernel()
    f=np.ascontiguousarray(f,dtype=np.float64);valid=np.ascontiguousarray(valid,dtype=np.uint8)
    out=np.zeros((len(f),18));components=np.full(len(f),-1,dtype=np.int32);active=components.copy()
    assert KERNEL.analyse(len(f),q,dim,f,valid,float(eps),float(top),float(floor),out,components,active)==0
    return out,components,active

def sample_walks(f,valid,eps,starts,greedy,seed,q=20,dim=4,steps=STEPS):
    starts=np.ascontiguousarray(starts,dtype=np.int32);paths=np.zeros((len(starts),steps+1),dtype=np.int32)
    assert KERNEL.walks(len(f),q,dim,np.ascontiguousarray(f),np.ascontiguousarray(valid,dtype=np.uint8),starts,len(starts),steps,float(eps),int(greedy),int(seed),paths)==0
    return paths

def validate_kernel(output=None):
    """Independent recursive enumeration and sparse components on small graphs."""
    checks=0
    for q,dim in [(2,3),(3,2)]:
        n=q**dim;codes=np.array(np.unravel_index(np.arange(n),(q,)*dim)).T
        neighbours=[np.flatnonzero((codes!=codes[u]).sum(axis=1)==1) for u in range(n)]
        for trial in range(5):
            rng=np.random.default_rng(6800+trial);f=np.round(rng.normal(size=n),1);valid=np.ones(n,dtype=bool)
            if trial%2:valid[1]=False
            for eps in [0.,.15]:
                target=float(np.quantile(f[valid],.8));floor=0.;out,comp,ac=graph(f,valid,eps,target,floor,q,dim)
                ranks=np.zeros(n);ranks[valid]=(rankdata(f[valid],method='average')-.5)/valid.sum();best=max(f[valid])
                @functools.lru_cache(None)
                def solve(u,greedy):
                    vv=[v for v in neighbours[u] if valid[v] and f[v]-f[u]>eps]
                    if greedy and vv:vv=[v for v in vv if f[v]==max(f[w] for w in vv)]
                    if not vv:return np.array([0,ranks[u],float(f[u]>=target),float(f[u]==best),f[u]])
                    z=np.mean([solve(v,greedy) for v in vv],axis=0);z[0]+=1
                    if f[u]>=target:z[2]=1
                    if f[u]==best:z[3]=1
                    return z
                for u in np.flatnonzero(valid):
                    for greedy,cols in [(False,[3,5,7,9,16]),(True,[4,6,8,10,17])]:np.testing.assert_allclose(out[u,cols],solve(int(u),greedy),atol=1e-12)
                    np.testing.assert_equal(out[u,0],sum(valid[neighbours[u]]))
                    for tol,col in [(eps,11),(eps,12),(None,13),(None,14)]:
                        targets=set(np.flatnonzero(valid & (f>=target if col in [11,13] else f==best)));seen={int(u)};queue=[int(u)]
                        for x in queue:
                            for v in neighbours[x]:
                                if valid[v] and v not in seen and ((f[v]>=f[x]) if tol is None else (f[v]-f[x]>tol)):seen.add(int(v));queue.append(int(v))
                        assert out[u,col]==bool(seen&targets)
                for selection,cc in [(valid,comp),(valid&(f>floor),ac)]:
                    ix=np.flatnonzero(selection)
                    adj=np.array([[int(v in neighbours[u] and abs(f[v]-f[u])<=eps) for v in ix] for u in ix])
                    if len(ix):
                        _,labels=connected_components(csr_matrix(adj),directed=False)
                        np.testing.assert_array_equal(labels[:,None]==labels[None,:],cc[ix,None]==cc[ix][None,:])
                checks+=1
    result={'small_graph_cases':checks,'independent_checks':['recursive walk moments','exhaustive reachability','SciPy components','degrees'],'kernel_sha256':sha(HERE/'landscape_graph_kernel.cpp')}

    if output is not None:dump(Path(output)/'kernel_validation.json',result)
    return result

def component_summary(f,components,valid,floor):
    selected=(components>=0);ids,counts=np.unique(components[selected],return_counts=True)
    if not len(ids):return {'n':0,'largest_fraction':None,'largest_span':None,'largest_below_floor_fraction':None,'largest_size':0,'mean_size_seen_by_node':None}
    biggest=ids[np.argmax(counts)];values=f[components==biggest]
    return {'n':len(counts),'largest_fraction':float(counts.max()/len(f[valid])),'largest_span':float(np.ptp(values)),
      'largest_below_floor_fraction':float((values<=floor).mean()),'largest_size':int(counts.max()),'mean_size_seen_by_node':float((counts.astype(float)**2).sum()/counts.sum())}
