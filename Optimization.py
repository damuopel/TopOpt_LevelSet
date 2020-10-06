# Classes
class TopOpt():
    def __init__(self,x,xold,v,r,p):
        self.x = x
        self.xold = xold
        self.v = v
        self.r = r
        self.p = p
        self.c = 0.0
        self.dc = np.zeros(x.shape)
        self.dv = np.ones(x.shape)
         
    def CreateFilter(self,Mesh):
        self.H = np.zeros((Mesh.nElms,Mesh.nElms))
        coords = Mesh.xy[:,Mesh.Topology.T.flatten()]
        x = np.sum(coords[0,:].reshape(Mesh.nElms,4),axis=1)/4
        y = np.sum(coords[1,:].reshape(Mesh.nElms,4),axis=1)/4
        center = np.array([x,y])
        for iElm in range(0,Mesh.nElms):
            for jElm in range(0,Mesh.nElms):
                d = np.sum(np.subtract(center[:,jElm],center[:,iElm])**2,axis=0)**0.5
                self.H[jElm,iElm] = np.maximum(0,self.r-d)
    
    def FilterSensitivities(self):
        sumH = np.sum(self.H,axis=0)
        self.dc = self.H@(self.x*self.dc)/(self.x*sumH[:,None])
    
    def OptimalityCriteria(self,Mesh):
        # Bisection Method
        l1 = np.amin(np.absolute(self.dc))/np.amax(np.absolute(self.dv))
        l2 = np.amax(np.absolute(self.dc))/np.amin(np.absolute(self.dv))
        move = 0.2
        while abs((l2-l1)/(l2+l1))>=tol:
            lmid = 0.5*(l1+l2)
            xnew = np.maximum(1e-3,np.maximum(self.x-move,np.minimum(1.0,np.minimum(self.x+move,self.x*np.sqrt(-self.dc/(self.dv*lmid))))))
            if np.sum(xnew) - self.v*Mesh.nElms > 0:
                l1 = lmid
            else:
                l2 = lmid
        self.x = xnew