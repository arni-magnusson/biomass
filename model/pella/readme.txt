2014, MAIN CHANGES

pella2011 - estimates logp
            p is between 0 and Inf
            Bmsy/K >= 0.368

pella2014 - estimates p
            user sets the bounds p in model.ctl, somewhere between -1 and Inf
            Bmsy/K can be anything, only restricted by p bounds in model.ctl

The model.ctl file needs to be changed to set sensible starting value and bounds
for logp and p.

fletcher  - different parametrization: m and n instead of r and p



OTHER DIFFERENCES

More output: Harvest rate, surplus production, reference points, etc.
