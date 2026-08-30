#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Genera las figuras numéricas usadas en el cierre de la presentación.

Sistemas:
1) Oscilador Yang--Mills--Higgs homogéneo de dos modos
   H = 1/2(p1^2+p2^2) + a(q1^2+q2^2) + 1/2 q1^2 q2^2.
   a=1 reproduce la normalización gráfica usada para el sistema tipo Salasnich;
   a=1/8 representa el coeficiente g^2 v_EW^2/8 con g=v_EW=1 de la ruta de doblete.
2) Yang--Mills homogéneo puro: a=0, potencial x^2 y^2/2.
3) Sistema efectivo de Canfora--Grandi--Oyarzo--Oliva:
   V(G,W)=2W^2G^2 + W^4/2 + lambda/4 (G^2-nu^2)^2,
   arXiv:2309.04915, ecs. (2.20)--(2.23).

Las secciones de Poincaré se construyen registrando cruces con q2=0, p2>0
(o G=0,p_G>0 para Canfora et al.).
"""
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from numba import njit

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "images" / "figures"
OUT.mkdir(parents=True, exist_ok=True)

# Colores deliberadamente de alto contraste para proyección en sala iluminada.
NAVY = "#082F49"
BLUE = "#0B5FA5"
RED = "#B42318"
ORANGE = "#C45A00"
GREEN = "#146C43"
BLACK = "#111111"
GRAY = "#4B5563"

plt.rcParams.update({
    "font.size": 10,
    "axes.titlesize": 10,
    "axes.labelsize": 9,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "axes.edgecolor": BLACK,
    "axes.linewidth": 0.9,
    "figure.facecolor": "white",
    "axes.facecolor": "white",
    "savefig.facecolor": "white",
})

@njit(cache=True)
def grad_quartic(q1, q2, a):
    return 2.0*a*q1 + q1*q2*q2, 2.0*a*q2 + q2*q1*q1

@njit(cache=True)
def V_quartic(q1, q2, a):
    return a*(q1*q1 + q2*q2) + 0.5*q1*q1*q2*q2

@njit(cache=True)
def poincare_quartic(E, a, q1_init, p1_init, dt, nsteps, max_points):
    norb = q1_init.size
    q1 = q1_init.copy()
    q2 = np.zeros(norb)
    p1 = p1_init.copy()
    p2 = np.empty(norb)
    active = np.ones(norb, dtype=np.bool_)
    for j in range(norb):
        rem = 2.0*(E - V_quartic(q1[j], 0.0, a)) - p1[j]*p1[j]
        if rem <= 0.0:
            active[j] = False
            p2[j] = 0.0
        else:
            p2[j] = np.sqrt(rem)
    xs = np.empty((norb, max_points))
    ys = np.empty((norb, max_points))
    counts = np.zeros(norb, dtype=np.int64)
    for _ in range(nsteps):
        old_q1 = q1.copy(); old_q2 = q2.copy(); old_p1 = p1.copy(); old_p2 = p2.copy()
        for j in range(norb):
            if not active[j] or counts[j] >= max_points:
                continue
            g1,g2 = grad_quartic(q1[j],q2[j],a)
            ph1 = p1[j]-0.5*dt*g1
            ph2 = p2[j]-0.5*dt*g2
            nq1 = q1[j]+dt*ph1
            nq2 = q2[j]+dt*ph2
            ng1,ng2 = grad_quartic(nq1,nq2,a)
            np1 = ph1-0.5*dt*ng1
            np2 = ph2-0.5*dt*ng2
            q1[j]=nq1; q2[j]=nq2; p1[j]=np1; p2[j]=np2
            if old_q2[j] < 0.0 and nq2 >= 0.0 and np2 > 0.0:
                den = nq2-old_q2[j]
                theta = 0.0 if abs(den)<1e-15 else (-old_q2[j])/den
                k = counts[j]
                xs[j,k] = old_q1[j] + theta*(nq1-old_q1[j])
                ys[j,k] = old_p1[j] + theta*(np1-old_p1[j])
                counts[j] += 1
        done=True
        for j in range(norb):
            if active[j] and counts[j] < max_points:
                done=False; break
        if done: break
    return xs,ys,counts

def quartic_initials(E, a, norb=9, seed=2):
    rng=np.random.default_rng(seed)
    if a>0:
        qmax=np.sqrt(E/a)*0.88
    else:
        qmax=max(1.4,0.85*np.sqrt(2*E)+1.0)
    q1=np.linspace(-qmax,qmax,norb)
    # pequeñas variaciones de p1 para muestrear varias órbitas sobre la misma energía
    pscale=np.sqrt(2*E)
    p1=0.10*pscale*np.sin(np.linspace(0,2*np.pi,norb,endpoint=False))
    # rescala donde sea necesario para que la energía cinética restante sea positiva
    for j in range(norb):
        avail=2*(E-V_quartic(q1[j],0,a))
        if avail<=0:
            q1[j]*=.75; avail=2*(E-V_quartic(q1[j],0,a))
        lim=np.sqrt(max(avail,1e-12))*0.65
        p1[j]=np.clip(p1[j],-lim,lim)
    return q1.astype(np.float64),p1.astype(np.float64)

def plot_quartic_panel(ax,E,a,title,max_points=900,norb=9,dt=0.012,nsteps=420000):
    q1,p1=quartic_initials(E,a,norb=norb,seed=int(100*E+10*a))
    xs,ys,c=poincare_quartic(E,a,q1,p1,dt,nsteps,max_points)
    for j in range(len(c)):
        if c[j]>20:
            ax.scatter(xs[j,:c[j]],ys[j,:c[j]],s=2.2,c=NAVY,alpha=.88,linewidths=0)
    ax.set_title(title,fontweight='bold')
    ax.set_xlabel(r'$q_1$')
    ax.set_ylabel(r'$p_1$')
    ax.axhline(0,lw=.45,c="#777777",alpha=.55)
    ax.axvline(0,lw=.45,c="#777777",alpha=.55)
    return int(c.sum())

# ---------- Lyapunov: método de dos trayectorias con renormalización ----------
@njit(cache=True)
def rhs_quartic(x,a):
    q1,q2,p1,p2=x
    g1,g2=grad_quartic(q1,q2,a)
    return np.array((p1,p2,-g1,-g2))

@njit(cache=True)
def rk4_step(x,dt,a):
    k1=rhs_quartic(x,a)
    k2=rhs_quartic(x+0.5*dt*k1,a)
    k3=rhs_quartic(x+0.5*dt*k2,a)
    k4=rhs_quartic(x+dt*k3,a)
    return x+(dt/6.0)*(k1+2*k2+2*k3+k4)

@njit(cache=True)
def running_lyapunov(E,a,q1_0,p1_0,dt,nsteps,renorm_steps,eps):
    p2sq=2*(E-V_quartic(q1_0,0,a))-p1_0*p1_0
    p2=np.sqrt(max(p2sq,1e-14))
    x=np.array((q1_0,0.0,p1_0,p2))
    d=np.array((0.37,-0.21,0.51,0.74))
    d=d/np.sqrt((d*d).sum())*eps
    y=x+d
    nout=nsteps//renorm_steps
    ts=np.empty(nout); lam=np.empty(nout); sep=np.empty(nout)
    acc=0.0; k=0
    for i in range(nsteps):
        x=rk4_step(x,dt,a); y=rk4_step(y,dt,a)
        if (i+1)%renorm_steps==0:
            d=y-x; dist=np.sqrt((d*d).sum())
            if dist<1e-300: dist=1e-300
            acc += np.log(dist/eps)
            t=(i+1)*dt
            ts[k]=t; lam[k]=acc/t; sep[k]=dist
            d=d/dist*eps; y=x+d
            k+=1
    return ts,lam,sep

# ---------- Canfora et al. ----------
@njit(cache=True)
def V_canfora(G,W,lam,nu):
    return 2*W*W*G*G + 0.5*W**4 + 0.25*lam*(G*G-nu*nu)**2
@njit(cache=True)
def grad_canfora(G,W,lam,nu):
    return 4*W*W*G + lam*G*(G*G-nu*nu), 4*W*G*G + 2*W**3
@njit(cache=True)
def poincare_canfora(energies,lam,nu,pW0,dt,nsteps,max_points):
    n=energies.size
    G=np.zeros(n); W=np.zeros(n); pW=np.full(n,pW0); pG=np.empty(n)
    for j in range(n): pG[j]=np.sqrt(max(2*energies[j]-pW0*pW0,1e-14))
    xs=np.empty((n,max_points)); ys=np.empty((n,max_points)); counts=np.zeros(n,dtype=np.int64)
    for _ in range(nsteps):
        oldG=G.copy(); oldW=W.copy(); oldpG=pG.copy(); oldpW=pW.copy()
        for j in range(n):
            if counts[j]>=max_points: continue
            gG,gW=grad_canfora(G[j],W[j],lam,nu)
            phG=pG[j]-.5*dt*gG; phW=pW[j]-.5*dt*gW
            nG=G[j]+dt*phG; nW=W[j]+dt*phW
            ngG,ngW=grad_canfora(nG,nW,lam,nu)
            npG=phG-.5*dt*ngG; npW=phW-.5*dt*ngW
            G[j]=nG; W[j]=nW; pG[j]=npG; pW[j]=npW
            if oldG[j] < 0.0 and nG >= 0.0 and npG>0.0:
                den=nG-oldG[j]; theta=0.0 if abs(den)<1e-15 else (-oldG[j])/den
                k=counts[j]
                xs[j,k]=oldW[j]+theta*(nW-oldW[j])
                ys[j,k]=oldpW[j]+theta*(npW-oldpW[j])
                counts[j]+=1
        done=True
        for j in range(n):
            if counts[j]<max_points: done=False; break
        if done: break
    return xs,ys,counts

# ---- 1. Poincaré principal: seis energías ----
energies=[2.5,5.5,6.2,8.0,10.0,14.0]
fig,axs=plt.subplots(2,3,figsize=(11.2,6.3),constrained_layout=True)
for ax,E in zip(axs.flat,energies):
    n=plot_quartic_panel(ax,E,1.0,rf'$E={E:g}$',max_points=1050,norb=10)
fig.suptitle(r'Secciones de Poincaré: $q_2=0$, $p_2>0$',fontsize=13,fontweight='bold')
fig.savefig(OUT/'salasnich_poincare_6E_contraste.png',dpi=300,bbox_inches='tight')
plt.close(fig)

# Intro: una regular y una caótica, en grande
fig,axs=plt.subplots(1,2,figsize=(8.8,3.7),constrained_layout=True)
plot_quartic_panel(axs[0],2.5,1.0,r'Región regular: $E=2.5$',max_points=1200,norb=8)
plot_quartic_panel(axs[1],10.0,1.0,r'Región mixta/caótica: $E=10$',max_points=1200,norb=8)
fig.savefig(OUT/'poincare_intro_real.png',dpi=300,bbox_inches='tight')
plt.close(fig)

# 2. Exponente de Lyapunov: estimador corrido para dos energías
curves=[]
# Elegimos deliberadamente una órbita regular y una órbita caótica del mismo Hamiltoniano.
for E,q10,p10 in [(2.5,0.0,0.0),(10.0,-0.8,-0.8)]:
    ts,lamv,_=running_lyapunov(E,1.0,q10,p10,0.005,220000,40,1e-8)
    curves.append((E,ts,lamv))
fig,ax=plt.subplots(figsize=(7.5,4.1),constrained_layout=True)
for (E,t,l),col,label in zip(curves,[GREEN,RED],['órbita regular, $E=2.5$','órbita caótica, $E=10$']):
    ax.plot(t,l,lw=2.7,c=col,label=label)
ax.axhline(0,c=BLACK,lw=.9)
ax.set_xlabel('tiempo')
ax.set_ylabel(r'estimación de $\lambda_{\max}(t)$')
ax.set_title('Exponente de Lyapunov máximo: estimación numérica',fontweight='bold')
ax.legend(frameon=True,loc='upper right')
ax.set_xlim(0,curves[0][1][-1])
fig.savefig(OUT/'lyapunov_exponente.png',dpi=300,bbox_inches='tight')
plt.close(fig)

# 3. Dos reducciones adicionales: doblete EW y Yang--Mills puro
fig,axs=plt.subplots(2,2,figsize=(8.6,6.8),constrained_layout=True)
for ax,E in zip(axs[0],[0.5,3.0]):
    plot_quartic_panel(ax,E,0.125,rf'$E={E:g}$',max_points=900,norb=9,dt=.012,nsteps=420000)
for ax,E in zip(axs[1],[0.5,2.0]):
    plot_quartic_panel(ax,E,0.0,rf'$E={E:g}$',max_points=850,norb=8,dt=.008,nsteps=380000)
axs[0,0].text(-.18,1.25,'SU(2)-Higgs (doblete, $g=v_{EW}=1$)',transform=axs[0,0].transAxes,
              fontsize=11,fontweight='bold',color=BLUE)
axs[1,0].text(-.18,1.25,'Yang-Mills homogéneo puro ($a=0$)',transform=axs[1,0].transAxes,
              fontsize=11,fontweight='bold',color=ORANGE)
fig.savefig(OUT/'reducciones_adicionales_poincare.png',dpi=300,bbox_inches='tight')
plt.close(fig)

# Figura compacta del sector electrodébil homogéneo para la diapositiva comparativa.
fig,axs=plt.subplots(1,2,figsize=(7.8,3.5),constrained_layout=True)
for ax,E in zip(axs,[0.5,3.0]):
    plot_quartic_panel(ax,E,0.125,rf'$E={E:g}$',max_points=1100,norb=10,dt=.012,nsteps=450000)
fig.suptitle(r'Sector $SU(2)$--Higgs homogéneo, $g=v_{EW}=1$',fontsize=11,fontweight='bold')
fig.savefig(OUT/'electroweak_poincare_contraste.png',dpi=300,bbox_inches='tight')
plt.close(fig)

# 4. Canfora--Grandi--Oyarzo--Oliva, lambda=1, nu=0, pW(0)=0.3818.
# El paper reporta transición alrededor de E=0.127--0.128 para el caso nu=0.
ener=np.array([0.123,0.127,0.128,0.135],dtype=np.float64)
xs,ys,c=poincare_canfora(ener,1.0,0.0,0.3818,0.01,1500000,1800)
fig,axs=plt.subplots(2,2,figsize=(8.6,6.6),constrained_layout=True)
for j,(ax,E) in enumerate(zip(axs.flat,ener)):
    ax.scatter(xs[j,:c[j]],ys[j,:c[j]],s=2.6,c=NAVY,alpha=.9,linewidths=0)
    ax.set_title(rf'$E={E:.3f}$',fontweight='bold')
    ax.set_xlabel(r'$W$'); ax.set_ylabel(r'$p_W$')
    ax.axhline(0,lw=.45,c="#777777",alpha=.55); ax.axvline(0,lw=.45,c="#777777",alpha=.55)
fig.suptitle(r'Canfora--Grandi--Oyarzo--Oliva: sección $G=0$, $p_G>0$',fontsize=12,fontweight='bold')
fig.savefig(OUT/'canfora_poincare_contraste.png',dpi=300,bbox_inches='tight')
plt.close(fig)

# 5. Aumentar contraste de la figura de diones BPS ya disponible en el proyecto.
# No se inventan puntos nuevos: se oscurecen y engrosan los marcadores existentes para proyección.
try:
    from PIL import Image
    from scipy.ndimage import binary_dilation
    src=OUT/'dyons_poincare_original.png'
    if src.exists():
        im=np.asarray(Image.open(src).convert('RGB')).copy()
        r,gc,b=im[:,:,0],im[:,:,1],im[:,:,2]
        mask=(b>120)&(gc>70)&(r<100)&((b-r)>55)
        mask=binary_dilation(mask,iterations=1)
        im[mask]=np.array([8,47,73],dtype=np.uint8)
        proc=Image.fromarray(im)
        proc.save(OUT/'dyons_poincare_contraste.png')
        # Paneles representativos: cuasiperiódico (arriba izq.) y caótico (abajo izq.).
        w,h=proc.size
        top=proc.crop((0,0,w//2,h//2+35))
        bot=proc.crop((0,h//2-25,w//2,h))
        tw=max(top.width,bot.width); th=max(top.height,bot.height)
        compact=Image.new('RGB',(tw*2,th),'white')
        compact.paste(top,(0,0)); compact.paste(bot,(tw,0))
        compact.save(OUT/'dyons_poincare_2panel_contraste.png')
except Exception as exc:
    print('Aviso: no se pudo procesar diones:',exc)

print('Figuras generadas en',OUT)
