# Version 2026 (2026-04-02)

* Added surplus production calculation for the last year in report file.




# Version 2014 (2014-04-28)

* The p parameter is now estimated on a linear scale, rather than a log scale.
  This allows p to have values between -1 and Inf, so Bmsy/K can range from 0 to
  1, rather than being restricted between 0.368 and 1. The user sets the bounds
  for the p parameter in model.ctl.

* Added reference points, harvest rate, and surplus production to report file.

* Added fletcher.tpl alternative parmetrization, estimating m and n instead of r
  and p.




# Version 2011 (2011-01-06)

* Improved vector extraction.

* Improved string handling and comments.




# Version 2010 (2010-03-09)

* Initial version, for benchmarking against R optimizers.
