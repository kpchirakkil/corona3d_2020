# Corona 3D Monte‑Carlo Hot‑Atom Transport Model Updates

The Corona 3D code simulates **hot‐H** and **hot‐O** coronae and photochemical escape at Mars and Venus.  
This file lists the concrete action items required to (i) integrate the latest **doubly differential** and **inelastic** collision cross‑sections and (ii) run the updated model with revised MAVEN in‑situ data set.

## Purpose
This document tracks action items for updating the Monte‑Carlo hot‑O escape model so it reflects the latest physics (new doubly differential & inelastic cross‑sections) and the revised MAVEN in‑situ data set.

## Implementation including doubly differential collision cross-sections

- Make sure the energy and angular dependent O-CO2 collision cross-sections are incorporated into the model correctly.
- Agreed that **elastic‑only treatments underestimate energy loss**; incorporating state‑resolved inelastic channels via branching algorithm will provide physically realistic escape probabilities.  
- Priority cross‑sections are **O–CO₂** (elastic + inelastic), followed by **O–O, O–N₂, O–CO**.
- **Collision branching algorithm**
  1. For each collision energy *E*, decide if a collision occurs using **total** σ(E).  
  2. Select **collision channel** (elastic vs. inelastic₁,₂,…) weighted by their partial cross‑sections.  
  3. Sample outgoing energy & scattering angle from the appropriate doubly differential σ(E,θ).
- Add **energy‑bin lookup** (≈ 0.5 eV spacing) to accelerate σ queries.  
- Keep legacy “elastic‑only” mode behind a `--elastic-only` flag for benchmarking.  
- Implement unit tests covering energy conservation, detailed balance, and limiting cases.
- Cross‑Section Data
    - **Elastic (doubly differential)**: Ingest σ(E,θ) for **O–CO₂, O–O, O–N₂, O–CO** (Gacesa et al. (2020) & companion datasets).  
    - **Inelastic**: Import state‑resolved σ(E,θ) from the Gacesa GitHub repository
        - O-CO2: <https://github.com/mgacesa66/O-CO2_cross-sections>
        - O-CO: <https://github.com/mgacesa66/Cross-sections-O-CO>
        - O-N2: <https://github.com/snchtchhbr/n2_o_cross_section>

## Task Checklist and Recent Discussion

- **Cross‑Section Data**
  - Gather doubly differential *elastic* cross‑sections for **O‑CO₂, O‑O, O‑N₂, O‑CO** from Gacesa et al. (2020) & related databases.
  - Collect corresponding *inelastic* cross‑sections for the same pairs (from the Gacesa GitHub repository: O-CO2 at <https://github.com/mgacesa66/O-CO2_cross-sections>, O-CO at <https://github.com/mgacesa66/Cross-sections-O-CO>, and O-N2 at <https://github.com/snchtchhbr/n2_o_cross_section>).
  - Convert all cross‑sections into interpolation tables (energy × scattering‑angle).

- **Monte‑Carlo Code**
  - Implement branching algorithm: 
  1. decide if a collision occurs (total σ)
  2. choose process (elastic vs. inelastic_i)
  3. sample outgoing energy & angle via differential σ.
  - Introduce energy‑bin handling (e.g., 0.5 eV spacing) to speed σ look‑ups.
  - Preserve previous “elastic‑only” mode behind a feature flag for benchmarking.
  - Unit‑test new collision kernel with analytic limits (energy conservation, detailed balance).

- **Need for Inelastic Collisions**  
  Highlighted that treating inelastic events as elastic underestimates energy loss; proposed a collision‑branch algorithm to capture proper energy & angular redistribution.

- **Cross‑Section Priorities**  
  O‑CO₂ data (elastic & inelastic) from Gacesa et al. (2020) are highest priority; next are O‑O, O‑N₂, O‑CO.

- **Validation**
  Reproduce MAVEN Deep‑Dip cases with new physics; compare to original elastic‑only runs.

- **Comparison Strategy**  
  The upgraded 3‑D Monte‑Carlo will be benchmarked against Lillis et al. (2017) MAVEN in-situ analysis, using the revised MAVEN mission data set.

- **Automation Goal**  
  Automate MC runs for many MAVEN orbits to map eacape flux dependence on season, SZA and EUV with the new physics.

- **Documentation & Release**
  Update comments and `README.txt`.

- **References**  
  - Gacesa et al. 2020, *MNRAS* **491**, 5650 (elastic & inelastic O–CO₂ σ(E,θ)) <https://academic.oup.com/mnras/article/491/4/5650/5651174>
  - Lillis et al. 2017, *JGR* **122**, 3815 (MAVEN hot‑O escape analysis) <https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/2016JA023525>

## Notes

- Divide the energy into a few regions instead of treating it as a continuous variable with respect to determining the elastic vs inelastic scattering behavior. For example, at 2.5-3 eV, we could take one average cross section to represent it. There will be differential cross sections, elastic and inelastic, for it, and they could be used in the code. I am guessing that treating the energy change of the particle is not a bad idea. However, the largest energy changes (kinetic -> internal) occur for large scattering angles, which means a large change in the direction. 
One approach just treated inelastic collisions as though they were elastic, where one assumes that the escaping particle cannot gain much energy from excited background particles (because they de-excite fast enough? – but the environment is non-thermal…)
- Seems the best way to include inelastic collisions in the code would be to think of them as branches of a generic collision process:
1. For the collision energy, look up the total cross section to determine if a collision has occurred
2. Determine the type of collision (inelastic processes 1,2,3,... vs. elastic) using the relative size of the cross sections at that energy
3. Once the type of collision is known, determine the resulting particle energy and scattering from the angular / energy differential cross section of the elastic or inelastic process.
Depending on what kind of information we actually have for the cross sections we would need to make defensible simplifying assumptions for any of the above.

