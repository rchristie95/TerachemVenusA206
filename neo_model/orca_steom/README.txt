DLPNO-STEOM-CCSD in-protein excitation (free, memory-safe alternative to ADC(2)).
Everything is prepared EXCEPT the ORCA binary (free, but needs YOUR academic account).

WHAT'S HERE
  geom_cthrp.xyz  field.pc       : chromophore + Tyr203-phenol (44 atoms) -> compare to ADC(2) 451 nm
  geom_big.xyz    field_big.pc   : chromophore + His148 + FULL Tyr203 (73 atoms) -> production
  steom_{svp,tzvp}.inp           : small model, def2-SVP (method comparison) / def2-TZVP (quality)
  steom_big_{svp,tzvp}.inp       : big model
  run_orca.sh                    : runner

ONE STEP YOU MUST DO (I can't make an academic account):
  1. Register (free) at https://orcaforum.kofo.mpg.de  (academic email)
  2. Download "ORCA 6.x, Linux, x86-64, shared-version" + the matching OpenMPI it names
  3. tar -xf the ORCA archive to e.g. /home/robson/orca ; install that OpenMPI
THEN I (or you) run:
  ORCABIN=/home/robson/orca ./run_orca.sh steom_svp.inp     # ~comparison first
  ORCABIN=/home/robson/orca ./run_orca.sh steom_big_tzvp.inp # production

WHY THIS FITS WHERE ADC(2) DIDN'T: DLPNO local correlation -> tiny memory; 73 atoms/def2-TZVP
runs in <~40 GB. STEOM-CCSD >= ADC(2)/CC2 quality. Closed-shell singlet (OK). Point-charge field included.
