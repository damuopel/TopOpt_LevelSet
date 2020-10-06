# Classes
class FEM():        
    def K(self,Mesh,Material,Top):
        # Integration
        # Initialize some variables
        xyGPI = np.array([[-0.5774,-0.5774,0.5774,0.5774],[-0.5774,0.5774,-0.5774,0.5774]])
        hGPI = np.array([1,1,1,1])
        refNodes = Mesh.Topology[:,0] # Pick a reference element (are the same)
        refVerts = Mesh.xy[:,refNodes]
        Ke = 0
        for iGP in range(4):
            Xi = xyGPI[0,iGP]
            Eta = xyGPI[1,iGP]
            H = hGPI[iGP]
            dNl = ShapeFunctions(Xi,Eta,1)
            Jacobian = refVerts@dNl.T
            dNg = inv(Jacobian)@dNl
            B = np.array([[dNg[0,0],0,dNg[0,1],0,dNg[0,2],0,dNg[0,3],0],[0,dNg[1,0],0,dNg[1,1],0,dNg[1,2],0,dNg[1,3]],[dNg[1,0],dNg[0,0],dNg[1,1],dNg[0,1],dNg[1,2],dNg[0,2],dNg[1,3],dNg[0,3]]])
            BtD = np.dot(B.T,Material.D)
            Ke = Ke + t*np.dot(BtD,B)*det(Jacobian)*H  
        # Assembly
        elsDofs = 2*np.kron(Mesh.Topology,np.ones((2,1)))+np.tile(np.array([[0],[1]]),[4,Mesh.nElms])
        row = np.kron(elsDofs,np.ones((1,dofsElements))).T.flatten()
        col = np.kron(elsDofs,np.ones((dofsElements,1))).T.flatten()
        K0 = np.tile(Ke.flatten(),(Mesh.nElms,1))
        data = np.repeat(Top.x,dofsElements**2,axis=1)**Top.p*K0
        K = csc_matrix((data.flatten(),(row,col)),shape=(Mesh.dofs,Mesh.dofs))
        return K,K0
    
    def Loads(self,Mesh):
        self.F = np.zeros((Mesh.dofs,1))
        if NuemannCase == 'Puntual Force':
            # Punctual Force
            maxValues = Mesh.xy.max(1)
            xMax = Mesh.xy[0,:]==maxValues[0]
            minValues = Mesh.xy.min(1)
            yMin = Mesh.xy[1,:]==minValues[1]
            nodes = np.arange(Mesh.nNodes)
            forceNodes = nodes[np.where(np.logical_and(xMax,yMin))]
            forceDofs = np.array([2*forceNodes,2*forceNodes+1])
            fx = 0
            fy = -1
            self.F[forceDofs[0,:]] = fx
            self.F[forceDofs[1,:]] = fy
            
    def Solver(self,Mesh,K):
        totalDofs = np.arange(Mesh.dofs)
        if DirichletCase == 'Clamp':
            minValues = Mesh.xy.min(1)
            xMin = Mesh.xy[0,:]==minValues[0]
            restNodes = np.where(xMin)[0]
            restDofs = np.array([2*restNodes,2*restNodes+1]).T.flatten()
            uRest = np.zeros((restDofs.size,1))
            freeDofs = np.setdiff1d(totalDofs,restDofs)
            
        freeDofsX,freeDofsY = np.meshgrid(freeDofs,freeDofs)
        frDofsX,frDofsY = np.meshgrid(restDofs,freeDofs)
        self.u = np.zeros((Mesh.dofs,1))
        uFree = linalg.inv(K[freeDofsX,freeDofsY])*self.F[freeDofs]-K[frDofsX,frDofsY]*uRest
        self.u[freeDofs] = uFree
        self.u[restDofs] = uRest
        return self.u

class Mesh():
    def __init__(self,nx,ny):
        self.nx = nx
        self.ny = ny
        self.nElms = nx*ny
        self.nNodes = (nx+1)*(ny+1)
        self.dofs = self.nNodes*dofsNode
        # Coordinate 
        x = np.linspace(0,h*nx,nx+1)
        y = np.linspace(0,h*ny,ny+1)
        x, y = np.meshgrid(x,y)
        self.xy = np.array([x.flatten(),y.flatten()])    
        # Topology
        n1 = np.array([iElm+floor(iElm/nx) for iElm in range(nx*ny)])
        n2 = n1 + 1
        n3 = n2 + nx + 1
        n4 = n1 + nx + 1
        self.Topology = np.array([n1,n2,n3,n4])    

class Material():
    def __init__(self):
        # Plain Stress
        self.D = (E/(1-nu**2))*np.array([[1,nu,0],[nu,1,0],[0,0,(1-nu)/2]])


# Functions
def ShapeFunctions(Xi,Eta,dFlag):
    if dFlag==1:
        dN1dXi = -0.25*(1-Eta)
        dN1dEta = -0.25*(1-Xi)
        dN2dXi = 0.25*(1-Eta)
        dN2dEta = -0.25*(1+Xi)
        dN3dXi = 0.25*(1+Eta)
        dN3dEta = 0.25*(1+Xi)
        dN4dXi = -0.25*(1+Eta)
        dN4dEta = 0.25*(1-Xi)
        N = np.array([[dN1dXi,dN2dXi,dN3dXi,dN4dXi],[dN1dEta,dN2dEta,dN3dEta,dN4dEta]])
    else:
        N1 = 0.25*(1-Xi)*(1-Eta)
        N2 = 0.25*(1+Xi)*(1-Eta)
        N3 = 0.25*(1+Xi)*(1+Eta)
        N4 = 0.25*(1-Xi)*(1+Eta)
        N = np.array([N1,N2,N3,N4])       
    return N