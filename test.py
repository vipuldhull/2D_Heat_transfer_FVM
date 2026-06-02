shape =  'quadratic'    # 'rectangular', 'linear',  'quadratic', 'crazy'

l = 1
dimX = 21
dimY = 21

boundary =   ['D', 'N', 'D', 'D'] # [N,S,W,E] : D : Dirichlet, N : Neumann, R : Robin
TD =  [100 , 100 , 300 , 100] # [N,S,W,E]
alpha = 20
Tinf = 90
q = 0

X, Y = setUpMesh(dimX, dimY, l, formfunction, shape)
heat = SteadyHeat2D_FVM(X, Y, boundary, TD, q, alpha, Tinf)
T = heat.solve()
heat.plot_Result(T, "2D")
heat.plot_Result(T, "3D")

#test case 2: Mixed BC
shape =  'linear'    # 'rectangular', 'linear',  'quadratic', 'crazy'

l = 1
dimX = 41
dimY = 41

boundary =   ['R', 'N', 'D', 'R'] # [N,S,W,E] : D : Dirichlet, N : Neumann, R : Robin
TD =  [100 , 100 , 100 , 100] # [N,S,W,E]
alpha = 20
Tinf = 90
q = 0

X, Y = setUpMesh(dimX, dimY, l, formfunction, shape)
heat = SteadyHeat2D_FVM(X, Y, boundary, TD, q, alpha, Tinf)
T = heat.solve()
heat.plot_Result(T, "2D")
heat.plot_Result(T, "3D")
