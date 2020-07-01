# Shelf Enabled Global NEMO (SE-NEMO)

**_\*\* NB the code as it stands is a placeholder - not intended for use - the setup script has been tested and will checkout and compile the code on ARCHER, but namelists and forcing files are yet to be configured_**

Configuration files for SE-NEMO project

## Quick Start:

```
git clone git@github.com:NOC-MSM/SE-NEMO.git
./SE-NEMO/scripts/setup/se-orca025_setup_archer -w $PWD/test -x $PWD/test -s $PWD/SE-NEMO
cd test/nemo/cfgs/se-orca025/EXP00
```
Edit the project code in  `runscript.pbs` then:
```
qsub -q short runscript.pbs
```

forcing data:

(could be held here)
[SE-ORCA025](http://gws-access.ceda.ac.uk/public/jmmp_collab/)

