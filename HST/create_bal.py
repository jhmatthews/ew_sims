dbal = np.empty([4,len(bal_indices)], dtype="string")

dbal[0] = BALs
d_hst["ra"][bal_indices]
d_hst["dec"][bal_indices]
d_hst["z"][bal_indices]

f = open("bal_catalog.dat", "w")

for i, b in enumerate(bal_indices):

	line = "%s %.5f %.5f %.5f\n" % (BALs[i], d_hst["ra"][b], d_hst["dec"][b], d_hst["z"][b]) 

	f.write(line)

f.close()