from ew_stats_main import *
from scipy.optimize import minimize

THMIN = 90.0
THMAX = 90.0

# get the data
data = sub.get_sample("../data/catalog.dat")

# selection criteria
# we need to be at a certain redshift and have O3 emission
# then we select for BAL or no-BAL
select = selection(data)


size_fit = int(1e6)

# data
ews = data["ew_o3"][select.general*select.nonbal]

costhetas, flags = get_mock_angles(THMIN, len(ews), THMAX)
#costhetas2, flags2 = get_mock_angles(THMIN, size_fit, THMAX)

min_obj = minimize(function_to_minimise, [10.0,5.0], args = (ews, costhetas),
                   bounds = ((1,50),(1,50)), method="Powell")

mu, sig = min_obj.x
chi2 = min_obj.fun
#mu1=mu2 = 8.1
#sig1 = 13.5
#sig2 = 6.

print mu, sig, chi2

ews_fake = np.random.normal(loc=mu, scale=sig, size=len(ews)) / costhetas
#ews_fake2 = np.random.normal(loc=mu2, scale=sig2, size=len(ews)) / costhetas

# make plot
bins = np.arange(0,200,1)
hist(ews, bins=bins, alpha = 0.4)
hist(ews_fake, bins=bins, alpha = 0.4)
savefig("hist.png")
#hist(ews_fake2, bins=bins, alpha = 0.4)

