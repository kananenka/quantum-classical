#
# Use this python script to generate .ff file
#
nwaters     = 500
outfilename = "water.ff"
h_charge    = 0.41
o_charge    = -0.82
h_mass      = 1.008
o_mass      = 15.9994
h_eps       = 0.0
o_eps       = 0.65
h_sig       = 0.0
o_sig       = 3.166

#
# prepare *.ff file
#
f = open(outfilename, "w")
f.write(" %d \n"%(3*nwaters))
for n in xrange(nwaters):
   f.write(" %d  \t %f \t  %f \t  %f \t %f \n"%(3*n,o_mass,o_charge,o_sig,o_eps))
   f.write(" %d  \t %f \t  %f \t  %f \t %f \n"%(3*n+1,h_mass,h_charge,h_sig,h_eps))
   f.write(" %d  \t %f \t  %f \t  %f \t %f \n"%(3*n+2,h_mass,h_charge,h_sig,h_eps))
f.close()

print " Force field file: ", outfilename, " has been generated."
