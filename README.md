# Aircraft Flight Dynamics Simulator (Fortran Physics + Unreal Engine Visualization)

Flight-dynamics simulator built for Aug–Dec 2025 (USU). Physics engine in modern Fortran with real-time visualization in Unreal Engine. Unreal receives state telemetry via UDP and renders aircraft motion/HUD in real time.

## Highlights
- Fortran physics engine: atmospheric model, nonlinear stall model, RK4 integrator, trimming
- Real-time Unreal visualization driven by physics state via UDP
- Config-driven runs (example inputs included)

## Tech Stack
**Fortran**, **Unreal Engine**, **UDP**

## Repo Structure
- `fortran/` – physics engine source + build + example inputs
- `unreal/` – Unreal setup notes + Blueprint/UDP integration notes
- `scripts/` – plotting / post-processing helpers
- `docs/` – methodology, equations, assumptions
- `media/` – simulator demo


