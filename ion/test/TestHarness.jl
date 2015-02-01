module TestHarness


Q = [ [ -3050 50           0      3000   0; ]
             [ 2/3   -(500 + 2/3) 500    0      0; ]
             [ 0     15000        -19000 4000   0; ]
             [ 15    0            50     -2065  2000;]
             [ 0     0            0      10     -10;] ]


tres = 2.5e-5
tcrit = 1e-4
conc = 1e-7
k=5
nopen=2


export Q,tres,tcrit,conc,k,nopen

end


