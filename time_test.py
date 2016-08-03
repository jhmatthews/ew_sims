import time
import ew_stats_main as ew
import numpy as np
import ew_sub

print ew_sub.__doc__

print "reading in..."
data = ew.sub.get_sample("data/catalog.dat")
sel = ew.selection(data)




fluxes = data["L5100"][sel.R11]
NPTS = 14424

print "all set up"

thmin = 5
thmax = 10
t1 = time.time()

#costhetas_qsos, dummy = ew.get_mock_angles2(0.0, NPTS, fluxes, max_angle=thmin)
costhetas1, mock_bal_flags1 = ew.get_mock_angles(thmin, NPTS, fluxes, max_angle=thmax)

t2 = time.time()

print t2 - t1

t1 = time.time()

#costhetas_qsos, dummy = ew.get_mock_angles(0.0, NPTS, fluxes, max_angle=thmin)
costhetas2, mock_bal_flags2 = ew.get_mock_angles2(thmin, NPTS, fluxes, max_angle=thmax)

t2 = time.time()

print t2 - t1

print "means:"
print np.mean(costhetas1), np.mean(costhetas2)