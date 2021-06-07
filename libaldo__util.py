from sympy import *   
from libaldo_polyclass  import *
from libaldo_lib import *
import lib_algebraEq as ae
 
numvar=-1
k0,x,a,b =symbols('k0 x a b')

vece=[k0,k0,k0,k0,k0,k0,k0,k0,k0,k0,k0,k0,k0]
e0,e1,e2,e3,e4,e5,e6,e7,e8,e9,e19,e11,e12=k0,k0,k0,k0,k0,k0,k0,k0,k0,k0,k0,k0,k0


   
 
    


def sqrs(k1,k2): # return pow(k1,S('k2'))
    v1=str(k2)
    v2='1/'+v1
    kres=pow(k1,S(v2))
    return(kres)
def frs(k1,k2): # return   S('k1/k2')
    v1=str(k1)
    v2=str(k2) 
    v3=v1+'/'+v2
    return(S(v3))

def rpow(ksym,op1='',op2='',kope=''):
    if type(op2)==int:
        ke=frs(op1,op2)
        kope=kope
    else:
        ke=frs(1,op1)
        kope=op2
    kres=(pow(ksym,ke))
    kres=opemat(kres,kope)
    return(kres)
    
def kpow(ksym,op1='',op2='',kope=''):
    if type(op2)==int:
        ke=frs(op1,op2)
        kope=kope
    else:
        ke=op1
        kope=op2
    kres=(pow(ksym,ke))
    kres=opemat(kres,kope)
    return(kres)

 
def get_num_mono(ksym,kveri=''):
    try:
        vlist=fpoly(ksym,'free')
        
        kres=ksym
        for i in vlist:
            kres=kres.subs(i,1)
        if kres==1 and kveri!='':
            return(False)
        elif kres!=1 and kveri!='':
            return(True)
        else:
            return(kres) 
    except:
        if kveri!='':
            return(True)
        else:    
         return(ksym)
        
        
def get_var_mono(ksym,kveri=''):
    try:
        vlist=fpoly(ksym,'free')
        kres1=ksym
        for i in vlist:
            kres1=kres1.subs(i,1)
        kres=ksym/kres1    
        if kveri!='':
            return(True)
        else:
            return(kres)
    except:
        if kveri!='':
            return(False)
        else:
            return(1)

   
 

    

    
def perf_root(ksym,op1=''): # try reducc
    x1,y1=symbols('x1,y1', rational=True,positive=True)
    kres=ksym
    if op1=='32':
        x1=monoroot(fpoly(kres,'get',0))
        y1=monoroot(fpoly(kres,'get',1))
        if ksym==expand(pow(x1+y1,2)):
            kres=(pow(x1+y1,2))
        elif ksym==expand(pow(x1-y1,2)):
            kres=(pow(x1-y1,2))    
    return(kres)

def monoroot(ksym):
    kk=ksym
    try:
        ww=fpoly(ksym,'list')
        mm=[]
        for i in ww:
            if type(i)==Pow:
                w1= fpoly(i,'get',1)
                w2=opemat(w1/2,'s')
                w3=fpoly(i,'get',0)
                w4=pow(w3,w2)
                mm.append(w4)
            else:
                mm.append(sqrs(i,2))
        kk=1
        for i in mm:
            kk=kk*i
    except:
        kk=sqrs(ksym,2)
    return(kk)    


def fpoly(ksym,kopt='',op2='',op3=''):
    '''
     'n': return number of args   
     'list': return list of args
     'get': return args No =op2 in list args
     'get_inv': return inverse of args No=op2 in,list
     'gets': return sum of term(a)+term(b)+... op2= 'ab...'
     'getp': return multi  of term(a)+term(b)+... op2= 'ab...'
     'free': return list of all variable symbols in ksym
     'filt':  
     'unfilt':  
     'frac2sum':  
     'lam':  return lambdify(op2,ksym)
     'simpow': 
     'subs': return(ksym.subs(op2,op3)
     'facmono': factorize ksym whit op2
     'in_allexp':
     'simp_fac':
     'simp_facs':
    '''
    karg=ksym.args
    klist=list(karg)
    knum=len(klist)
    kres=ksym
    done=False
    if kopt=='h':
        kres=kres.subs(op2,op3)
        print(op2,op3,kres.subs(op2,op3))
         
        
    if kopt=='n':
        kres=(knum)
        done=True
       
    if kopt=='list':
    
        kres=(klist)
        done=True
    if kopt=='get':
        kres=(klist[op2])
        done=True
    if kopt=='get_inv':
        kres=(klist[op2])
        kres=kres**-1
        done=True    
    if kopt=='gets':
        done=True
        for i in op2:
            nsym=klist[int(i)]
            if done:
                mm=nsym
                done=False
            else:
                mm=mm+nsym
        kres=(mm)        
        done=True 
    if kopt=='getp':
        mm=1
        for i in op2:
            nsym=klist[int(i)]
            mm=mm*nsym
        kres=(mm)        
        done=True     
    if kopt=='free':
         kres=(list(ksym.free_symbols ))
         done=True
    if kopt=='filt' and op2!='':
        kres=0
        for i in fpoly(ksym,'list'):
            if op2 in fpoly(i,'free'):
                kres=kres+i
    
    if kopt=='filtp' and op2!='':
        kres=0
        for i in fpoly(ksym,'list'):
            if op2 in fpoly(i,'free'):
                kres=kres+i    
                
    if kopt=='unfilt' and op2!='':
        kres=0
        for i in fpoly(ksym,'list'):
            if op2 not in fpoly(i,'free'):
                kres=kres+i             
         
         
    if kopt=='frac2sum': 
        kres=klist
         
        kq0=klist[0]
        kq1=kq0.args[0]
        kq2=klist[1]
        kres=(kq2,kq1)
         
        done=True
    if kopt=='rqt':
        kstr='1/'+str(op2)
        kres=Pow(ksym,S(kstr))
    if kopt=='lam':
        kres=lambdify(op2,ksym)
    if kopt=='simpow': 
        kres=pow(fpoly(ksym,'get',0).args[0],fpoly(ksym,'list')[0].args[1]*fpoly(ksym,'get',1))    
    if kopt=='zubs':
        kres=ksym.subs(op2,op3)
    if kopt=='zubsV':
        for i, j in zip(op2,op3):
            kres=kres.subs(i,j)
             
    if kopt=='facmono':
        klist=poly(ksym,'list')
        mm=0
        for i in klist:
            if op2 in fpoly(i,'free'):
                mm=mm+i
        kres=mm    
    if kopt=='in_allexp':
        bres=True
        klist=fpoly(ksym,'list')
        for i in klist:
            if op2 not in fpoly(i,'free'):
                bres=False
        kres=bres
    
    if kopt=='simp_fac':
        newm=0
        oldm=0
        klist=fpoly(ksym,'list')
        kvvar=fpoly(ksym,'free')
        kvar=fpoly(op2,'free')
        kvar=kvar[0]
        for i in klist:
            kres1=fpoly(i/op2,'free')
            if kvar not in kres1:
                newm=newm+i/op2
            else:
                oldm=oldm+i
        if op3==1:
            kres=(newm)
        elif op3==2:
            kres=(oldm)
        elif op3==0:
            kres=(op2*newm)
        else:    
            kres=op2*(newm)+oldm
    
    if kopt=='simp_facs':
         
        kres=ksym
        skres=0
        veck=op2
        for op2 in veck:

            klist=fpoly(kres,'list')
            km1=0
            km2=0
            for i in klist:
                try: 
                    mm= fpoly(i,'list')
                    if len(mm)>0:
                        if op2 in mm:
                            #print(mm)
                            km1=km1+i/op2
                            km2=km2+i
                except:
                    done=False
            kres=kres-km2
            skres=skres+op2*km1
        skres=skres+(kres) 
        kres=skres
        
    if kopt=='list_tree':
        mm=[]
        kres=mm
        for i in fpoly(ksym,'list'):
            mm.append(short_type_name(i))                                  #mm.append(short_type_name(i))
        kres=mm     
        
    if kopt=='get_type':
         kres=short_type_name(ksym)
        
       
           
    return(kres)
 
def simple_fac(ksym,vecv):
    kres=ksym
    sexp=0
    for i in vecv:
        exp1=fpoly(kres,'facmono',i) 
        exp2=fpoly(exp1,'simp_fac',i) 
        sexp=sexp+exp2
        kres=kres-exp1 
    return(sexp)    
     
       
    
    
def polyargs(ksym,kargs):
    kres=ksym
    sexp=0
    for i in kargs:
        numar=int(i)
        kres=kres.args[numar]
    return(kres)

    
def pao(ksym,kargs):
    return(polyargs(ksym,kargs)) 


def paolist(ksym,kop=''):
    done=True
    cc=0
    while done:
        try:
            newb=kop+str(cc)
            kres=pao(ksym,newb)
            show_res(kres,str(cc))
            cc+=1
        except:
            done=False
        
    
def argolist(ksym,kv=''): 
    try:    
        mm=fpoly(ksym,'list')
        nn=0
        for i in mm:
            nkv=kv+str(nn)
            get_args_list(i,nkv)
            nn+=1
            show_res(i,nkv)  
    except:
        show_res(ksym,kv)      
      
# TRANSFORMACIONE CON lambdify


def polytrans(op1='',k1='',k2='',k3='',kdis=''):
    x,y,z=symbols('x y z')
    
    if op1=='h':
        
        show_res(ae.Equation((x+y)*(x**2-x*y+y**2),(x**3+y**3)),'op9')
        show_res(ae.Equation((x**3+y**3),(x+y)*(x**2-x*y+y**2)),'op10')
        show_res(ae.Equation(x**3+y**3+z**3-3*x*y*z,(x+y+z)*(x**2+y**2+z**2-x*y-y*z-z*x)),'op11')
        
    else:
        x=k1
        y=k2
        z=k3
        if op1=='9':
            kres=(x**3+y**3)
        if op1=='10':
            kres=(x+y)*(x**2-x*y+y**2) 
        if op1=='11':
            kres=(x+y+z)*(x**2+y**2+z**2-x*y-y*z-z*x) 

        return(kres)

def perfect_sqrt_exp(ksym,kop=''):
    klist=fpoly(ksym,'list')
    korden=([0,1,2],[0,2,1],[1,2,0])
    for kk in korden:
        e1=ksymroot(klist[kk[0]]) 
        e2=klist[kk[1]] 
        e3=ksymroot(klist[kk[2]])
    
        if 2*e1*e3==e2:
            kres=e1+e3
            if '1' in kop:
                kres=(e3+e1)
            if '2' in kop:
                kres=kpow(kres,2)
            return(kres)    
        elif -2*e1*e3==e2:
            kres=e1-e3
            if '1' in kop:
                kres=(e3-e1)
            if '2' in kop:
                kres=kpow(kres,2)
            return(kres)    
                
        else:
            kres=ksym 
    return(kres)        

def ksymroot(ksym): #perfect_sqrt_exp
    kres=1
    mm=gek(ksym)
    for i in mm:
        kv=i[0]
        ke=i[1]
        kres=kres*kpow(kv,ke,2)
    return(kres)

def gek(ksym): # used by ksymroot
    veriv=get_var_mono(ksym,1)
    verin=get_num_mono(ksym,1)
    mm=[]
    if verin:
        knum=get_num_mono(ksym)
        mm.append([knum,1])
    if veriv:
        kvar=get_var_mono(ksym)
        if type(kvar)==Symbol:
            mm.append([kvar,1])
        else:    
            if type(kvar)==Pow:
                mm.append([kvar.args[0],kvar.args[1]])
            else:
                klist=fpoly(kvar,'list')
                for i in klist:
                    if type(i)==Pow:
                        mm.append([i.args[0],i.args[1]])
                    else:
                        mm.append([i,1])
                    
    return(mm)                
    
    
def simpify_monoexp(ksym):
    try:  # metodo 1
        mm=fpoly(ksym,'list')
        if len(mm)==2:
            m1=mm[0]
            e1=mm[1]
            m1=opemat(m1,'s')
            mm1=m1.args[0]
            mm2=m1.args[1]
            mm3=opemat(e1*mm2,'s')
            return(kpow(mm1,mm3))
    except:
        return(ksym)
        
        
def simp_exp(ksym=1,op1='h'):
    if op1=='h':
        show_exp(['op=1,',kpow(kpow(x,a),b),'=',kpow(x,a*b)])
        show_exp(['op=2,',kpow(x,a)*kpow(x,b),'=',kpow(x,a+b)])
    elif op1==1:
            if type(ksym)==Pow:
                try:
                    mm=fpoly(ksym,'list')
                    ke1=mm[1] 
                    e2=mm[0] 
                    ke2=e2.args[1] 
                    m1=e2.args[0] 
                    kres=kpow(m1,opemat(ke1*ke2,'s'))
                    return(kres)
                except:
                    return(ksym)
            else:
                return(ksym)
    elif op1==2:
        if type(ksym)==Mul:
            try:
                mm=fpoly(ksym,'list')
                v1=mm[0] 
                v2=mm[1] 
                vv1=v1.args[0] 
                ve1=v1.args[1] 
                vv2=v2.args[0] 
                ve2=v2.args[1] 
                ke=opemat(ve1+ve2,'s') 
                kres=kpow(vv1,ke)
                return(kres)
            except:
                return(ksym)
        else:
            return(ksym)
    else:        
        return(ksym)
    
def short_tree_exp(ksym): # 'STE(ksym)' short called function 
    kres=srepr(ksym)
     
    nomL=['Symbol','Integer','Add','Mul','Pow','Rational','Float','.','1','2','3','4','5','6','7','8','9','00','000','0000','00000','000000','0000000','00000000','000000000']
    nomC=['S','I','A','M','P','R','F','','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0']
    for i,j in zip(nomL,nomC):
        kres=kres.replace(i,j)
    kres=kres.replace(' ','')
    nomL=['00','000','0000','00000','000000','0000000','00000000','000000000']
    nomC=['0','0','0','0','0','0','0','0']
    for i,j in zip(nomL,nomC):
        kres=kres.replace(i,j)
    kres=kres.replace('00','0')
    kres=kres.replace(',precision=0','')
    kres=kres.replace('(0)','')
    kres=kres.replace("'","")
    return(kres)
 
def short_type_name(ksym): # 'TT(ksym)' short called function 
    Nsym=type(ksym)
    if Nsym == Symbol :
        return('S')
    elif Nsym == Integer :
        return('I')
    elif Nsym == Add :
        return('A')
    elif Nsym == Mul :
        return('M')
    elif Nsym == Pow :
        return('P')
    elif Nsym == Rational :
        return('R')
    elif Nsym == Float :
        return('F')
    else:
        return('N')    
        
def TT(ksym):
    return(short_type_name(ksym))
def STE(ksym):    
    return(short_tree_exp(ksym))
    
    

# Operaciones con Monomios    
def monomath(ksym='',op1='',kope=''):
    x,a,b,c=symbols('x a b c')
    kres=ksym
    if op1=='h' or ksym=='h' or ksym=='':
        show_exp(['all letters includinn in op1 will be ejecuted'])
        show_exp(['monomath(f(x),SeOf...) or Eq.polymath(SeOf..) '])
        show_exp(['O','=',kpow(kpow(x,a),b),'to',kpow(x,a*b)])
        show_exp(['S','=',kpow(x,a),'to',kpow(x,(b+c))])
        show_exp(['F','=',kpow(x,a),'to',kpow(x,(b*c))])
        show_exp(['e expand, f factor,  s implify, 1 only left Eq, 2 only right Eq'])
    else:
        for i in op1:
        
            if i=='O':
                try:

                    if type(kres)==Pow:
     
                        kres2,exp1=kres.args[0],kres.args[1]
                        try:
                            if type(kres2)==Pow:
                                
                                 
                                kres3,exp2=kres2.args[0],kres2.args[1]
                                newexp=exp1*exp2
                                newexp=opemat(newexp,kope)
                                kres=kpow(kres3,newexp)
                                 
                                 
                            
                        except:
                            done=False
                except:
                    done=False        
            if i=='S':
                try: 
                    if type(kres)==Pow:
                         
                        kres2,exp1=kres.args[0],kres.args[1]
                        if type(exp1)==Mul:
                             
                            exp2=expand(exp1)
                             
                            if type(exp2==Add):
                                kres=(kpow(kres2,exp2))
                    else:
                        done=False  
                except:
                    done=False                
            if i=='F':
                try :
                    kres2,exp1=kres.args[0],kres.args[1]
                    exp2=factor(exp1)
                    kres=(kpow(kres2,exp2))
                except:
                    done=False
                            
            kres=opemat(kres,i) 
        
        return(kres) 
    
def get_exp_in_mono(ksym,kvar):
    klist=fpoly(ksym,'list')
    for i in klist:
        klist2=fpoly(i,'list')
        if kvar in klist2:
            return(klist2[1])
    return(ksym)
    


       