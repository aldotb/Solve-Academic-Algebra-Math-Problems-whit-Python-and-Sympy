from sympy import *   
from libaldo_lib import *
from libaldo__util import *
import lib_algebraEq as ae
a,b,c,d =symbols('a b c d')
x,x1,x2,y,y1,y2,z,z1,z2,t=symbols('x x1 x2 y y1 y2 z z1 z2 t')
kp= symbols('kp')
 
class monoclass:
    def __init__(self,q=kp):
        self.Q=q
        self.V=[]
        self.H=kp
    
    
    def showq(self):
        return(self.Q)
    def re_eQ(self,ksym):
        self.update(ksym)
        kres=self.Q
        return(kres)
        
    def opemat(self,kope='',kd='0' ):
        kres= opemat(self.Q,kope) 
        if kd==1:
                self.update(kres) 
        return(kres)
        
    def update(self,ksym):
        self.H=self.Q
        self.Q=ksym
        
    def ope4( self,op1='',op2=0,kd='0',kope=''):
        kres=self.Q
        if op1 in 'SMDP':
            if op1=='S':
                kres=kres+op2
            if op1=='M':
                kres=kres*op2     
            if op1=='D':
                kres=kres/op2
            if op1=='P':
                kres=pow(kres,op2)
            try:
                kres=opemat(kres,kd)
            except :
                kres=opemat(kres,kope)
                 
            if kd==1:
                self.update(kres)        
                 
        return(kres)
    
    def get_exp(self):
        kres=self.Q
        if type(kres)==Pow:
            kres=kres.args[1]   
        return(kres)
        
    def subs(self,op1,op2,kd=''):
        kres=self.Q 
        kres=kres.subs(op1,op2)
        return(kres)        
    
class polyclass:
    def __init__(self,p1=kp,p2=kp):
        self.q1=p1
        self.q2=p2
        self.Q=ae.Equation(p1,p2)
        self.V=[]
        self.H=ae.Equation(p1,p2)
    
    
    def showq(self):
        return(self.Q)
    
    
    
    def re_eQ(self,p1,p2):
        kres=ae.Equation(p1,p2)
        self.update(kres)
        return(kres)
    
    def update(self,newQ):
        self.H=self.Q
        self.Q=newQ
        self.q1=newQ.args[0]
        self.q2=newQ.args[1]
        
    def opemat(self,kope='',op1='',op2=''):
        if op1=='1':
            kres=ae.Equation(opemat(self.q1,kope),self.q2)
        elif op1=='2':
            kres=ae.Equation(self.q1,opemat(self.q2,kope))
        else:
            kres=self.Q
            kres=opemat(kres,kope)
        if op1==1 or op2==1:
            self.update(kres)
        return(kres)
    
        
    def ope4( self,op1='',op2=0,kd='0',kope=''):
        kres=self.Q
        if op1 in 'SMDP':
            if op1=='S':
                kres=kres+op2
            if op1=='M':
                kres=kres*op2     
            if op1=='D':
                kres=kres/op2
            
            if op1=='L':
                kres=kres/op2
            
            if op1=='P':
                kres=pow(kres,op2)
            try:
                kres=opemat(kres,kd)
            except :
                kres=opemat(kres,kope)
                 
            if kd==1:
                self.update(kres)        
                 
        return(kres)
        
    def ope1( self,op1='',op2=0,kd='0',kope=''):
        kres=self.q1
        if op1 in 'SMDP':
            if op1=='S':
                kres=kres+op2
            if op1=='M':
                kres=kres*op2     
            if op1=='D':
                kres=kres/op2
            if op1=='P':
                kres=pow(kres,op2)
            try:
                kres=opemat(kres,kd)
            except :
                kres=opemat(kres,kope)
                 
            if kd==1:
                newEq=ae.Equation(kres,self.q2)
                self.update(newEq) 
                return(self.Q)                
                #return(kres)
        return(ae.Equation(kres,self.q2))
                
    def ope2( self,op1='',op2=0,kd='0',kope=''):
        kres=self.q2
        if op1 in 'SMDP':
            if op1=='S':
                kres=kres+op2
            if op1=='M':
                kres=kres*op2     
            if op1=='D':
                kres=kres/op2
            if op1=='P':
                kres=pow(kres,op2)
            try:
                kres=opemat(kres,kd)
            except :
                kres=opemat(kres,kope)
                 
            if kd==1:
                newEq=ae.Equation(self.q1,kres)
                self.update(newEq) 
                return(self.Q)                
                #return(kres)
        return(ae.Equation(self.q1,kres))
        
    def swap(self):
        p2=self.q1
        p1=self.q2
        kres=ae.Equation(p1,p2)
        self.Q=kres
        return(kres)
    
    def solve(self,op1):
        kres=self.q1-self.q2
        kres=csolve(kres,op1)
        return(kres)
    def get_exp(self,op1,kope=''):
        
        kres=self.Q
        kres=kres.args[int(op1)-1]
        if type(kres)==Pow:
            kres=kres.args[1]
        kres=opemat(kres,kope)    
        return(kres)
    
    def simpPow(self,op1):
        kres1=self.q1
        kres2=symbols(self.q2)
        if op1==1:
            kres1=poly(kres1,'simpow')
        if op1==2:
            kres=type(kres2)  
        #kres=ae.Equation(kres1,kres2) 
        return(kres)    
    
    def subs(self,op1,op2='',kd=''):
        kres=self.Q
        if op2=='' or op2=='up':
            kres1=(self.q1).subs(op1)
            kres2=(self.q2).subs(op1)
            kres=ae.Equation(kres1,kres2)
            if op2=='up':
                kd=1
             
        else:
            kres=kres.subs(op1,op2)

        if kd!='':
                self.update(kres)
        return(kres)        
        
    def subsubs(self,op1,op2,kd=''):
        kres=self.Q
        for k1,k2 in zip(op1,op2):
            kres=kres.subs(k1,k2)
        if kd!='':
                self.update(kres)
        return(kres)
        
    def evalue_if(self,op1,op2,kope=''):
        kres1=self.q1
        kres1=kres1.subs(op1,op2)
        kres1=opemat(kres1,kope)
        
        kres2=self.q2
        kres2=kres2.subs(op1,op2)
        kres2=opemat(kres2,kope)
        
        kres=(kres1==kres2)
        print(kres)
        return(ae.Equation(kres1,kres2))
        
    def simpFac(self,kd=''):
        m1,m2=self.q1,self.q2
        kfac= gcd(m1,m2)
        kres1=simplify(m1/kfac)
        kres2=simplify(m2/kfac)
        kres=ae.Equation(kres1,kres2)
        if kd!='':
                self.update(kres)
        return(kres)
    
    def simp_log_LE(self,kd=''):
        kres=ae.Equation(elel(self.q1),elel(self.q2))
        if kd!='':
                self.update(kres)
        return(kres)
        
    def simp_log_S2M(self,kd=''):
        kres=ae.Equation(opelog2(self.q1,'f'),opelog2(self.q2,'f'))
        if kd!='':
                self.update(kres)
        return(kres)
        
        
        
        
    