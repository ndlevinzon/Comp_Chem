#Based on [https://pubs.acs.org/doi/10.1021/ct400862k]

parm AACCGA_1.D2O.top
trajin md.rst

#Dihedral Force Constant * 1.0
scaledihedralk AACCGA_1.D20.top 1.0 1:194
parmwrite out AACCGA_1.D2O.H10.top
trajout AACCGA_1.D2O.H10.rst7 restart

#Dihedral Force Constant * 0.9
scaledihedralk AACCGA_1.D20.top 0.9 1:194
parmwrite out AACCGA_1.D2O.H09.top
trajout AACCGA_1.D2O.H09.rst7 restart

#Dihedral Force Constant * 0.8
parm AACCGA_1.D2O.top
scaledihedralk AACCGA_1.D20.top 0.8 1:194
parmwrite out AACCGA_1.D2O.H08.top
trajout AACCGA_1.D2O.H08.rst7 restart

#Dihedral Force Constant * 0.7
parm AACCGA_1.D2O.top
scaledihedralk AACCGA_1.D20.top 0.7 1:194
parmwrite out AACCGA_1.D2O.H07.top
trajout AACCGA_1.D2O.H07.rst7 restart

#Dihedral Force Constant * 0.6
parm AACCGA_1.D2O.top
scaledihedralk AACCGA_1.D20.top 0.6 1:194
parmwrite out AACCGA_1.D2O.H06.top
trajout AACCGA_1.D2O.H06.rst7 restart

#Dihedral Force Constant * 0.5
parm AACCGA_1.D2O.top
scaledihedralk AACCGA_1.D20.top 0.5 1:194
parmwrite out AACCGA_1.D2O.H05.top
trajout AACCGA_1.D2O.H05.rst7 restart

#Dihedral Force Constant * 0.4
parm AACCGA_1.D2O.top
scaledihedralk AACCGA_1.D20.top 0.4 1:194
parmwrite out AACCGA_1.D2O.H04.top
trajout AACCGA_1.D2O.H04.rst7 restart

#Dihedral Force Constant * 0.3
parm AACCGA_1.D2O.top
scaledihedralk AACCGA_1.D20.top 0.3 1:194
parmwrite out AACCGA_1.D2O.H03.top
trajout AACCGA_1.D2O.H03.rst7 restart
