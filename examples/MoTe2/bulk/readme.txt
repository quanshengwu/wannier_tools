An example to obtain band structure of MoTe2 with the Hamiltonian obtained from the
OpenMX with non-orthogonal basis. 

The Hamiltonain is stored as two parts in sparse format, H.dat and S.dat.

H.dat contains the Hamiltonian Hmn(R)=<\phi_m(r)|H|\phi_n(r-R)>, only the non-zero entries are stored.
S.dat contains the overlap matrix Smn(R)=<\phi_m(r)|\phi_n(r-R)>, only the non-zero entries are stored

tar xzvf HS.tar.gz 
