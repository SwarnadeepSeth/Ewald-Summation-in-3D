import numpy as np

def coulomb_force(r1, r2, q1, q2):
	r_vec = r1-r2
	rsq = np.linalg.norm(r_vec)
	return q1*q2/(rsq**2)
	 
	 
fi = open("bare_coulomb.dat", "w")
z=1.0
while (z < 35):
	r1 = np.array([0.0,0.0,0.0])
	r2 = np.array([0.0,0.0,z])
	
	q1=1.0
	q2=-1.0
	
	f = coulomb_force(r1, r2, q1, q2)
	f1 = coulomb_force(r1, np.array([0.0,0.0,1.0]), q1, q2)
	
	print (z, f/f1, file = fi)
	
	z = z + 0.1

fi.close()
