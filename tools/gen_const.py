#
# Use this python script to generate .const file
#
nwaters     = 500
outfilename = "water.const"
rOH         = 1.0
ang_HOH     = 109.47

#
# prepare *.const file
#
n_bond  = 2*nwaters
n_angle = nwaters
f = open(outfilename, "w")
f.write(" %d %d \n"%(n_bond, n_angle))
for n in xrange(nwaters):
   f.write(" %d %d %f \n"%(3*n, 3*n+1, rOH))
   f.write(" %d %d %f \n"%(3*n, 3*n+2, rOH))
for n in xrange(nwaters):
   f.write(" %d %d %d %f \n"%(3*n, 3*n+1, 3*n+2, ang_HOH))
f.close()

print " Constraints file: ", outfilename, " has been generated."
