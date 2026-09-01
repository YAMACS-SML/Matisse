# Matisse
This repository contains all the supplimentary informations of Matisse
MATISSE, a novel algorithm that analyses protein binding pockets using water molecules and short molecular dynamics (MD) simulations. This method addresses the limitations of conventional tools by exploring the binding pocket's spatial and volumetric properties. MATISSE enhances our understanding of protein-ligand interactions by accounting for the dynamic geometry and flexibility of binding pockets.
Supporting Information

Item
Value / definition
Source in macro
Input structure
complex containing protein as object 1 and ligand as object 2
Loadcomplex; object checks
Measurement cycles
5
numeroCicli = 5
Simulation cell
Automatic cuboid; 1 Å extension
Cell Auto, Extension=1, Shape=Cuboid
Water density
0.997 g cm⁻³
FillCellWater Density=0.997
Solvent probe
1.4 Å
Probe=1.4
Force field
AMBER03
ForceField AMBER03
Nonbonded cutoff
10.5 Å
Cutoff 10.50000
Long-range treatment
None in the macro
Longrange None
Cavity-water cutoff
r < 4.1 Å to nearest protein atom
if raggioSfera < 4.1
Nonbonded force scaling during cycle
1.5
ScaleForce NonBonded,1.5
Pseudocode
ALGORITHM S1  MATISSE volumetric-topology characterization
 
INPUT:
protein–ligand structure with protein and ligand separately identifiable
N_cycles = 5
R_cut = 4.1 Å
OUTPUT:
V, stV, V_intersect, stV_intersect
 
A. STRUCTURE PREPARATION
1. Load the protein–ligand structure.
2. Remove residues named DOD.
3. Verify that a ligand object is present and contains more than five atoms.
4. Perform the structure-cleaning operation used for the production calculation.
5. Create an automatic cuboid simulation cell extending 1 Å beyond the system.
 
B. FIRST HYDRATION AND SOLVENT RELAXATION
6. Fill the cell with explicit water at density 0.997 g cm^-3 using a 1.4 Å probe.
7. Hold the protein and ligand objects fixed.
8. Use AMBER03 parameters, a 10.5 Å cutoff, and no long-range correction.
9. Enable bond, angle, dihedral, planarity, Coulomb, and van der Waals interactions.
10. Apply the macro's steepest-descent solvent-relaxation stage.
 
C. LIGAND REMOVAL AND SECOND HYDRATION
11. Remove the ligand from the hydrated complex.
12. Refill the same simulation cell with water using the same density and probe settings.
13. Identify the water population introduced by this second fill.
14. Let W = {O_1, O_2, ..., O_N} be the oxygen atoms of this second-fill water population.
 
D. INITIAL WATER–WATER DISTANCE MATRIX (AS IMPLEMENTED IN THE MACRO)
15. For every ordered pair (O_i, O_j) in W:
D_OO[i,j] = EuclideanDistance(O_i, O_j)
16. Store D_OO before the five measurement cycles.
 
E. FIVE ANNEALING / MEASUREMENT CYCLES
17. FOR cycle = 1 TO N_cycles:
17.1 Set temperature control to annealing.
17.2 Scale nonbonded forces by 1.5.
17.3 Continue the simulation for the waiting interval defined by the production setup.
17.4 Pause the simulation before geometric measurements.
17.5 Set V_cycle = 0 and retained_spheres = empty.
 
17.6 FOR each oxygen O_i in W:
P_i = nearest protein atom to O_i
r_i = EuclideanDistance(O_i, P_i)
IF r_i < R_cut:
v_i = (4/3) * pi * r_i^3
store sphere S_i = (water identity O_i, center O_i, radius r_i, volume v_i)
V_cycle = V_cycle + v_i
 
17.7 Set V_intersect_cycle_raw = 0.
17.8 FOR every ordered pair of stored radius entries (a,b):
r_a = stored radius for entry a
r_b = stored radius for entry b
d_ab = D_OO[a,b]   # pre-cycle matrix in the original macro
IF d_ab < (r_a + r_b) AND d_ab != 0:
V_ab = pi*(r_a+r_b-d_ab)^2 *
[d_ab^2 + 2*d_ab*r_a + 2*d_ab*r_b - 3*r_a^2 - 3*r_b^2 + 6*r_a*r_b]
/ (12*d_ab)
V_intersect_cycle_raw = V_intersect_cycle_raw + V_ab
17.9 V_intersect_cycle = V_intersect_cycle_raw / 2
# division by 2 removes double counting of (a,b) and (b,a)
17.10 Store V_cycle and V_intersect_cycle.
 
F. FINAL DESCRIPTORS
18. V = mean(V_cycle over the five cycles)
19. stV = standard_deviation(V_cycle over the five cycles)
20. V_intersect = mean(V_intersect_cycle over the five cycles)
21. stV_intersect = standard_deviation(V_intersect_cycle over the five cycles)
22. Return V, stV, V_intersect, stV_intersect.
Mathematical definitions used by the macro
Sphere radius: r_i = min_j ||O_i - P_j||
Individual sphere volume: v_i = (4/3)πr_i³
Cycle total sphere volume: V_cycle = Σ_i v_i
Two-sphere intersection volume: V_ij = π(r_i+r_j-d_ij)²[d_ij²+2d_ijr_i+2d_ijr_j-3r_i²-3r_j²+6r_ir_j] / (12d_ij)
Cycle pairwise-intersection sum: V_intersect,cycle = Σ_{i<j} V_ij
Downstream descriptor used in the manuscript. The uploaded macro returns V, stV, V_intersect, and stV_intersect. If the manuscript defines the effective-volume descriptor as V_e = V - V_intersect, this calculation is downstream of the macro and should be documented as such.
Non-proprietary Python reimplementation: suggested component mapping
The pseudocode is intentionally software-independent. A researcher could implement the same algorithmic operations using freely available Python packages. 
