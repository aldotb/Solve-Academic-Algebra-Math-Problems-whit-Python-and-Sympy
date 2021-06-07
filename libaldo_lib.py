from sympy import *
from IPython.display import display, Math
from libaldo_polyclass  import *
from libaldo_lib import *
import lib_algebraEq as ae

def multishow(vecs):
    for eeq in vecs:
        show_eq(eeq)
        
def show_eq(Eq1):
    kres=Eq1
    display(Math(latex(kres)))
    
def show_text(kd=''):
    sR=kd
    display(Math(sR))

def show_res(Eq1,kd='',kr=''):
    kres=Eq1
    sR=kd+' ='
    display(Math(sR+latex(kres)))
    if kr!='':
        return(kres)
        
def show_if(Eq1):
    kres=Eq1
    sR='if  \; \; \; \;' 
    display(Math(sR+latex(kres)))
def show_ifthen(Eq1):
    kres=Eq1
    sR1='if  \; \; \; \;'
    sR2='\; \; \; \; then'   
    display(Math(sR1+latex(kres)+sR2))
       
def fraseconvert(lst):
    return (lst[0].split())
def show_expP(vecexp):
    kexp=''
    for i in vecexp:
        if type(i)==str:
            kexp=kexp+ i + '\;'
        else:
            kexp=kexp+latex(i)+'\;'
    display(Math(kexp))
    
def show_exp(vecexp):
    kexp=''
    for i in vecexp:
        if type(i)==str:
            vecfrase=fraseconvert([i])
            for j in vecfrase:
                kexp=kexp+j + '\;'
        else:
            kexp=kexp+latex(i)+'\;'
    display(Math(kexp))     

def define_eQ(kstr,ksym):
    kexp=[]
    for vs,vv in zip(kstr,ksym):
        kexp.append(vs+'=')
        kexp.append(vv)
        kexp.append(', ')
        kexp.append('\; \; \;')
        kexp.append(' ')
    show_exp(kexp)
    return(ksym)    
        
        
        
def show_eques(vecvec):
    kexp=''
    for i in vecvec:
        kexp=str(i)+' \;'+'='+' \;'+latex(i)+', \;'
    display(Math(kexp))
    
def Ksym(Eq1,kd=''):
    kres=Eq1
    sR=kd+'='
    if kd!='':
        display(Math(sR+latex(kres)))
    return(kres)

def eQ(Eq1,kd=''):
    sR=kd+' ='
    display(Math(sR+latex(Eq1)))
    return(Eq1)        
 
def sex2rad(k):
    k=simplify(pi*k/180)
    return(k)

def rad2sex(k):
    k=simplify(180*k/pi)
    return(k.evalf())

def sex2rad_i(kang,s='r'):
    if s=='s':
        kang=sex2rad(kang)
    return(kang)
    

def csolve(Eq,ksym,kd='',korden='',kpositive=False,kope='',kdisp=False):
    if not(kd==''):
        kdisp=True
    
            
    kres=solve(Eq,ksym)
    qq=len(kres)
    if qq==1:
        kres=kres[0]
    if qq>1:
        kres=list(kres)
    if kpositive:
        kres=[i for i in kres if i>0]
        
        if len(kres)==1:
            kres=kres[0]
    if not(korden==''):
        if qq>1:
            kres=kres[korden]
            
    kres=opemat(kres,kope)
    #qq=len(kres)
    #if qq==1:
        #kres=kres[0]
    if kd!='':
        sR=kd+' ='
        display(Math(sR+latex(kres)))        
    #if kreturn:    
    return(kres)
        
        
def csolveset(Eq,kv,kd='',kdisp=False):
    kres=solveset(Eq,kv)
    kres=list(kres)
    if kdisp or kd!='':
        sR=kd+' ='
        display(Math(sR+latex(kres)))
    return(kres)    


        
 
       
def opemat(kres,kope=''):
    '''
    opemat(Equation,kope=opt)
    opt 
        'f'= factor(Equation),    
        't'= trigsimp(Equation)
        'x'= simplify(Equation)
        'v'= Equation.evalf()
        'a'= apart(Equation)
        'c'= cancel(Equation)
        'E'= Equation.expand(force=True)
    '''
    
    
    if 'e' in kope:
        kres=expand(kres)
    if 'f' in kope:
        kres=factor(kres)    
    if 't' in kope:
        kres=trigsimp(kres)
    if 'x' in kope:
        kres=expand_trig(kres)    
    if 's' in kope:
        kres=simplify(kres)
    if 'v' in kope:
        
        if count_symbols(kres)==0:
            try:
                kres=kres.evalf()
            except:  
                try: 
                    if not(type(kres)==float or type(kres)==int):
                        kres=kres.evalf()
                except:
                    kres=kres
                
                    
                    
    if 'a' in kope:
        kres=apart(kres)
    if 'c' in kope:
        kres=cancel(kres)
    if '2' in kope:
        kres=powsimp(kres**2)
        kres=powsimp(sqrt(kres))
    if 'E' in kope:
        kres=kres.expand(force=True)
    
    
    return(kres)        
    
def opemat_mat(kres,kope=''):
    mm=[]
    qq=len(kres)
    for i in range(qq):
        kkval=kres[i]
        kres1=opemat(kkval,kope)
        mm.append(kres1)
    return(mm)

def opemath(ksym,kope=''):
    kres=ksym
    for i in kope:
    
        if i=='D':
            kres=denom(kres)
        if i=='N':
            kres=numer(kres)	
        if i=='f':
            kres=factor(kres)    
        if i=='e':
            kres=expand(kres)
        if i=='s':
            kres=simplify(kres)
    
    return(kres) 
def ksimplify(kres,kope='',kd=' ',kdisp=False):
    if type(kres)==list:
        kres=opemat_mat(kres,kope)
    else:
        kres=opemat(kres,kope,kd,kdisp=kdisp)
        
    return(kres)
    
    
def clean_res(k,s='s'):
    kk=k[0]
    if s=='e':
        kk=float(kk)
    return(kk) 

def getv_ctable(kvar,ktable):
        kres=[]
        mm1=ktable
        qq=len(mm1)
        for i in range(qq):
            kv=mm1[i][0]
            if kv==kvar:
                kres=mm1[i]
        return(kres)    
def get_depen_func(kvar,ktable):
        kres=[]
        mm1=ktable
        qq=len(mm1)
        for i in range(qq):
            kv=mm1[i][0]
            if kv==kvar:
                kres=mm1[i][2]
        return(kres) 

def get_func(kvar,ktable,kdisp=False):
    kres=[]
    mm1=ktable
    qq=len(mm1)
    for i in range(qq):
        kv=mm1[i][0]
        if kv==kvar:
            kres=mm1[i][1]
    if kdisp:
        show_eq(kres)        
    return(kres) 

def string2symbol(a4):
    a5=[]
    qq=len(a4)
    qq
    for i in range(qq):
        kvec=a4[i]
        q2=len(kvec)
        a5v=[]
        for j in range(q2):
            kval=kvec[j]
            nkval=parse_expr(kval)
            a5v.append(nkval)
        a5.append(a5v)
    return(a5)

def get_matsec2var(ksym1,ksym2,matseq):
    qq=len(matseq)
    kres=[]
    for i in range(qq):
        kvec=matseq[i]
        if (ksym1 in kvec) and (ksym2 in kvec):
            kres.append(kvec)
    return kres

def diff2var(ksym1,ksym2,kmat):
    nfunc=get_func(ksym1,kmat)
    kres=nfunc.diff(ksym2)
    return(kres)

def diff_minivec(kvec,kmat,kdisp=False):
    mm2=kvec
    kk=1
    qq=len(mm2)
    for i in range(qq-1):
        var1=mm2[i]
        var2=mm2[i+1]
        kk=kk*diff2var(var1,var2,kmat)
    if kdisp:
        show_eq(kk)
    return(kk)

def diff_supervec(kvec,kmat,kdisp=False):
    mm2=kvec
    kres=0
    qq=len(mm2)
    for i in range(qq):
        kminv=mm2[i]
        kres=kres+diff_minivec(kminv,kmat,kdisp=False)
    return(kres)

def subssubs(ksym,vecvar,vecval,kmat):
    kres=[]
    kres=get_func(ksym,kmat)
    if kres==[]:
        kres=ksym
    qq=len(vecval)
    for i in range(qq):
        f1 = lambdify(vecvar[i], kres)
        kres=f1(vecval[i])
    return(kres) 
    
def get_args(keq,kk,kd=' ',kdisp=False):
    sR=kd+' ='
    kres=keq.args[kk]
    if kdisp or kd!=' ':
        display(Math(sR+latex(kres)))
    return(kres)

def get_sub2_arg(keq,k1,k2,kd=' ',kdisp=True):
    sR=kd+' ='
    kres=keq.args[k1]
    kres=kres.args[k2]
    if kdisp:
        display(Math(sR+latex(kres)))
    return(kres)

    
def count_args(keq):
    return( len(Eq.args))

    
def show_arg(keq,kfull=False):
    kres=[]
    mkarg=keq.args
    qq=len(mkarg)
    for i in range(qq):
        karg=mkarg[i]
        cc=i
        kres.append([i,karg])
    if kfull:    
        kres.append(['>>>',keq])
    display(Math(latex(kres)))
 
def get_subpoly_arg(keq,kk=[],kd=' ',kdisp=True):
    sR=kd+' ='
    qq=len(kk)
    kres=keq
    for i in range (qq):
        k1=kk[i]
        kres=kres.args[k1]
     
    if kdisp:
        display(Math(sR+latex(kres)))
    return(kres) 

def count_symbols(kvsym):
    kres=kvsym.free_symbols
    kres=list(kres)
    kqq=len(kres)
    return(kqq)
    
def simple_subs(kEq,ksym,kval,kope=''):
        kres=kEq.subs(ksym,kval)
        kres=opemat(kres,kope)
        return(kres)

def cdsolve(kEq,kfunc,ics='',kreturn=True):
    if ics!='':
        kres=dsolve(kEq,kfunc,ics=ics)
    else:
        kres=dsolve(kEq,kfunc,ics=ics)
    
    show_eq(kres)
    kres=kres.args[1]
    if kreturn:
        return(kres)
     
        
    
    
def simple_dsolveset(kEq,kfunc,kreturn=True):
    kres=dsolve(kEq,kfunc)
    kres1=kres[0]
    kres2=kres[1]
    show_eq(kres1)
    show_eq(kres2)
    if kreturn:
        return(kres[0].args[1],kres[1].args[1])
        
    

        
def sistem2obj(obj1,obj2):
    eq1=[obj1.ratio_eq(),obj2.ratio_eq()]
    eq2=[obj1.f,obj2.f]
    kres=dsolve(eq1,eq2)
    obj1.rsolu=kres[0].args[1]
    obj2.rsolu=kres[1].args[1]
    
# Geometric Functiones
def evald(self,kval,kope=''):
    kres=self.Dr
    kres=kres.subs(self.t,kval)
    kres=opemat(kres,kope)
    return(kres)   

def eq_tangent(self,x1,y1,kfull=False):
    
    eq1=self.y
    eq2=self.Df*(self.t-x1)+y1
     
    kres=Eq(eq1,eq2)
    if kfull:
        return(kres)
    else:    
        return(eq2)
    
def eq_normal(self,x1,y1,kfull=False): 
    eq1=self.f
    eq2=-1*(self.t-x1)/self.Df+y1
    kres=Eq(eq1,eq2)

    if kfull:
        return(kres)
    else:    
        kres=kres.args[1]
        return(eq2)

def long_tan(self,x1='',kope=''):
    if x1=='':
        x1=self.t
    kres=self.f*pow(1+(self.Dr**2),S('1/2'))/self.Dr
    kres=kres.subs(self.t,x1)
    kres=opemat(kres,kope)
    return(kres)

def long_subtan(self,x1='',kope=''):
    if x1=='':
        x1=self.t
    kres=self.f/self.Dr
    kres=kres.subs(self.t,x1)
    kres=opemat(kres,kope)
    return(kres)

def long_normal(self,x1='',kope=''):
    if x1=='':
        x1=self.t
    kres=self.f*pow(1+(self.Dr**2),S('1/2')) 
    kres=kres.subs(self.t,x1)
    kres=opemat(kres,kope)
    return(kres)

def long_subnorm(self,x1='',kope=''):
    if x1=='':
        x1=self.t
    kres=self.f*self.Dr 
    kres=kres.subs(self.t,x1)
    kres=opemat(kres,kope)
    return(kres)   

def despeja_diff(self,ksol):
    kres=csolve(ksol,self.Df)
    return(kres)

def logsolve(A='',B=0,C=0,D=0):
    if A=='':
        show_res(a*exp(b*k)+c,'d')
        show_text('logsolve(a,b,c,d)')
         
    else:
        kres=(ln((D-C)/A))/B
        return(kres)  
   