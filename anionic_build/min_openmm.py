from openmm.app import AmberPrmtopFile, AmberInpcrdFile, PME, HBonds, PDBFile, Simulation
from openmm import LangevinMiddleIntegrator, Platform
from openmm.unit import nanometer, kelvin, picosecond, picoseconds, kilojoule_per_mole
Platform.loadPluginsFromDirectory("/home/robson/anaconda3/envs/TeraChem/lib/plugins")
def pick():
    for n in ("CUDA","OpenCL","CPU","Reference"):
        try:
            p=Platform.getPlatformByName(n); return p, ({'Precision':'mixed'} if n in("CUDA","OpenCL") else {})
        except Exception: continue
prm=AmberPrmtopFile('monomer_solv.prmtop'); crd=AmberInpcrdFile('monomer_solv.inpcrd')
sys_=prm.createSystem(nonbondedMethod=PME, nonbondedCutoff=1.0*nanometer, constraints=HBonds)
integ=LangevinMiddleIntegrator(300*kelvin,1/picosecond,0.002*picoseconds)
plat,props=pick(); print("platform:",plat.getName(),flush=True)
sim=Simulation(prm.topology,sys_,integ,plat,props)
sim.context.setPositions(crd.positions)
if crd.boxVectors is not None: sim.context.setPeriodicBoxVectors(*crd.boxVectors)
e0=sim.context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(kilojoule_per_mole)
print(f"E(before) = {e0:.1f} kJ/mol",flush=True)
sim.minimizeEnergy(tolerance=10*kilojoule_per_mole/nanometer, maxIterations=20000)
st=sim.context.getState(getEnergy=True,getPositions=True)
print(f"E(after)  = {st.getPotentialEnergy().value_in_unit(kilojoule_per_mole):.1f} kJ/mol",flush=True)
PDBFile.writeFile(prm.topology, st.getPositions(), open('monomer_min.pdb','w'), keepIds=True)
print("wrote monomer_min.pdb",flush=True)
