Project-specific overlay environment guidance.

Rules:
1. Use the approved base stack by default.
2. Use requirements.in only for small project-specific additions.
3. Regenerate requirements.lock.txt after changing requirements.in.
4. Do not add GPU frameworks such as torch/JAX/CUDA here.
5. Do not commit virtual environments.
6. If a dependency is needed by multiple projects, request a registered stack update.
