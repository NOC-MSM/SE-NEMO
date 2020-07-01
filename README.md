# Shelf Enabled Global NEMO (SE-NEMO)

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

