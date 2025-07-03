# Corona 3D Monte‑Carlo Hot‑Atom Transport Model Updates

The Corona 3D code simulates **hot‐H** and **hot‐O** coronae and photochemical escape at Mars and Venus.  
This file lists the concrete action items required to (i) integrate the latest **doubly differential** and **inelastic** collision cross‑sections and (ii) re‑analyse the full MAVEN in‑situ data set.

## Purpose
This document tracks action items for updating the Monte‑Carlo hot‑O escape model so it reflects the latest physics (new doubly differential & inelastic cross‑sections) and the revised MAVEN in‑situ data set.

## Implementation including doubly differential collision cross-sections

- Make sure the energy and angular dependent O-CO2 collision cross-sections are incorporated into the model correctly.
- Agreed that **elastic‑only treatments underestimate energy loss**; incorporating state‑resolved inelastic channels via branching algorithm will provide physically realistic escape probabilities.  
- Priority cross‑sections are **O–CO₂** (elastic + inelastic), followed by **O–O, O–N₂, O–CO**.  A future effort will extend σ(E,θ) calculations to **N, C, H** collisions.
- **Collision branching algorithm**
  1. For each collision energy *E*, decide if a collision occurs using **total** σ(E).  
  2. Select **collision channel** (elastic vs. inelastic₁,₂,…) weighted by their partial cross‑sections.  
  3. Sample outgoing energy & scattering angle from the appropriate doubly differential σ(E,θ).
- Add **energy‑bin lookup** (≈ 0.5 eV spacing) to accelerate σ queries.  
- Keep legacy “elastic‑only” mode behind a `--elastic-only` flag for benchmarking.  
- Implement unit tests covering energy conservation, detailed balance, and limiting cases.
- Cross‑Section Data
    - **Elastic (doubly differential)**: Ingest σ(E,θ) for **O–CO₂, O–O, O–N₂, O–CO** (Gacesa 2020 & companion datasets).  
    - **Inelastic**: Import state‑resolved σ(E,θ) from the Gacesa GitHub repository <https://github.com/mgacesa66/O-CO2_cross-sections>.

##  Task Checklist and Recent Discussion

- **Cross‑Section Data**
  - Gather doubly differential *elastic* cross‑sections for **O‑CO₂, O‑O, O‑N₂, O‑CO** from Gacesa et al. (2020) & related databases.
  - Collect corresponding *inelastic* cross‑sections for the same pairs (from the Gacesa GitHub repository <https://github.com/mgacesa66/O-CO2_cross-sections>).
  - Convert all cross‑sections into interpolation tables (energy × scattering‑angle).

- **Monte‑Carlo Code**
  - Implement Mike’s branching algorithm: 
  1. decide if a collision occurs (total σ)
  2. choose process (elastic vs. inelastic_i)
  3. sample outgoing energy & angle via differential σ.
  - Introduce energy‑bin handling (e.g., 0.5 eV spacing) to speed σ look‑ups.
  - Preserve previous “elastic‑only” mode behind a feature flag for benchmarking.
  - Unit‑test new collision kernel with analytic limits (energy conservation, detailed balance).

- **Need for Inelastic Collisions**  
  Highlighted that treating inelastic events as elastic underestimates energy loss; proposed a collision‑branch algorithm to capture proper energy & angular redistribution.

- **Cross‑Section Priorities**  
  O‑CO₂ data (elastic & inelastic) from Gacesa et al. (2020): <https://academic.oup.com/mnras/article/491/4/5650/5651174> are highest priority; next are O‑O, O‑N₂, O‑CO. Future work will extend to N, C, H collisions.

- **Validation**
  Reproduce MAVEN Deep‑Dip cases with new physics; compare to original elastic‑only runs.

- **Comparison Strategy**  
  The upgraded 3‑D Monte‑Carlo will be benchmarked against Lillis et al. (2017): <https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/2016JA023525> MAVEN in-situ analysis, using the revised mission data set.

- **Automation Goal**  
  Automate MC runs for many MAVEN orbits to map eacape flux dependence on season, SZA and EUV with the new physics.

- **Documentation & Release**
  Update `README.md` and inline docstrings.

- **References**  
  - Gacesa et al. 2020, *MNRAS* **491**, 5650 (elastic & inelastic O–CO₂ σ(E,θ))  
  - Lillis et al. 2017, *JGR* **122**, 3815 (MAVEN hot‑O escape analysis)
