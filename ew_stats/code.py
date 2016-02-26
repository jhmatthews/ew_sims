from ew_stats_main import *
from scipy.optimize import minimize

THMIN = 20
THMAX = 20

distribution = np.random.lognormal

# get the data
data = sub.get_sample("../data/catalog.dat")

# selection criteria
# we need to be at a certain redshift and have O3 emission
# then we select for BAL or no-BAL
select = selection(data)


size_fit = int(1e6)

# data
ews = data["ew_c4"][select.b*select.nonbal]


costhetas, flags = get_mock_angles(THMIN, len(ews), THMAX)
#costhetas = np.ones(len(ews))
#costhetas2, flags2 = get_mock_angles(THMIN, size_fit, THMAX)

min_obj = minimize(function_to_minimize, [3.0,0.5], args = (ews, costhetas, distribution),
                   bounds = ((1,5),(1,5)), method="Powell")

mu, sig = min_obj.x
chi2 = min_obj.fun
#mu1=mu2 = 8.1
#sig1 = 13.5
#sig2 = 6.

print mu, sig, chi2

#ews_fake = np.random.normal(loc=mu, scale=sig, size=len(ews)) / costhetas
ews_fake = distribution(mu, sig, size=len(ews))
ews_fake = ews_fake / costhetas
#ews_fake2 = np.random.normal(loc=mu2, scale=sig2, size=len(ews)) / costhetas


bins = np.arange(0,200,1)

from pretty import *

set_pretty()
fig1=figure()
frame1=fig1.add_axes((.1,.3,.8,.6))
n1, bins1, patches1 = hist(ews, bins=bins, alpha = 0.4, label="EW$_O$")
n2, bins1, patches1 = hist(ews_fake, bins=bins, alpha = 0.4, label="EW$_{fake}$")

float_legend()
ylabel("$N$", fontsize=20)
chi = (n1 - n2) / np.sqrt(n1)
frame2 = fig1.add_axes((.1,.1,.8,.2))

print bins.shape, chi.shape
bins_plot = 0.5 * (bins[1:] + bins[:-1])
plot(bins_plot, chi, linewidth=2, c="k")
ylabel("$\chi$", fontsize=20)
xlabel("EW (\AA)", fontsize=20)
savefig("hist.png", dpi=200)
#hist(ews_fake2, bins=bins, alpha = 0.4)
clf()
close()

