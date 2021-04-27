"""
Utilities for working with BOUT++ output

Note, does not support all options in Python version:
* always includes all boundary cells
* variable names are case-sensitive
"""
module BoutPost

using BoutOutputs
export openboutdata
export boutcollect

using EndpointRanges
export ibegin, iend

end
