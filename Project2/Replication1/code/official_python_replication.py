"""Python port of the LIM simulations in deaneckles/eval_stochastic_seeding.

Sources (MIT-licensed project):
  simulations/main_simulations.R, village.R, outcomes.R, runner.R
This is an independently implemented Python port, not an author-provided script.
Results are statistically comparable, not random-seed identical to R.

Usage:
  python official_python_replication.py --raw data/0422allinforawnet.dta --replicates 10 --output outputs/official_quick
  python official_python_replication.py --raw data/0422allinforawnet.dta --replicates 5000 --output outputs/official_5000
"""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
from numba import njit
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ALPHAS=(-3.,-2.,-1.,0.)
BETAS=(0.,2.,4.,6.,8.,10.)
DEGREE_EFFECTS=(0.,0.1)


def load_author_villages(raw_file):
    """Port of main_simulations.R lines filtering villages and create_village_cai().

    The R code obtains vertices only from endpoints in graph_from_edgelist,
    filters on RAW edge row count > 25 and does not require survey data.
    Preserve parallel edges and self loops, as igraph graph_from_edgelist does.
    """
    raw=pd.read_stata(raw_file,convert_categoricals=False)
    good=raw.loc[(raw['address']==raw['network_address']) &
                 raw['network_id'].notna() & (raw['id']!=99) &
                 (raw['network_id']!=99) & raw['address'].notna()].copy()
    selected=good.groupby('address',sort=False).size()
    selected=selected[selected>25].index
    villages=[]
    diag=[]
    for addr in selected:
        subset=good.loc[good['address']==addr,['id','network_id']]
        nodes=pd.unique(subset[['id','network_id']].to_numpy().ravel())
        lookup={key:i for i,key in enumerate(nodes)}
        n=len(nodes)
        A=np.zeros((n,n),dtype=np.uint8)
        for a,b in subset.itertuples(index=False,name=None):
            A[lookup[a],lookup[b]]+=1
        outdeg=A.sum(axis=1,dtype=np.int64)
        valid=outdeg>0
        # R compute_p_nom_one_seed discards nominators with outdegree zero.
        q=(A[valid,:]/outdeg[valid,None]).mean(axis=0)
        assert abs(q.sum()-1)<1e-9
        norm=1-np.dot(q,q)
        if norm<=0:raise RuntimeError(f'No eligible seed pairs: {addr}')
        # R Village() uses degree(g, mode='in') for b_degree.
        indeg=A.sum(axis=0,dtype=np.int64)
        villages.append((str(addr),A,q,indeg,outdeg,norm))
        diag.append(dict(address=addr,nodes=n,raw_edges=len(subset),unique_edges=int(A.sum()),valid_nominators=int(valid.sum())))
    if len(villages)!=150:
        raise AssertionError(f'Expected 150 villages under author filters, got {len(villages)}. Check raw data version.')
    return villages,pd.DataFrame(diag)


@njit
def draw_one(q,n):
    r=np.random.random();c=0.
    for i in range(n):
        c+=q[i]
        if r<c:return i
    return n-1


@njit
def choose_seeds(q,n,onehop):
    if onehop:
        # Match R runner.R: choose an unordered pair with probability
        # 2*q[a]*q[b]/(1-sum(q**2)). Reject the ENTIRE ordered draw
        # when the two independent samples coincide.
        while True:
            a=draw_one(q,n)
            b=draw_one(q,n)
            if a!=b:
                break
    else:
        a=np.random.randint(n);b=np.random.randint(n-1)
        if b>=a:b+=1
    return a,b


@njit
def adoption(A,outdeg,indeg,n,seeda,seedb,alpha,b_degree,b_spill):
    y=np.zeros(n,np.uint8)
    y[seeda]=1;y[seedb]=1
    seed_degree=indeg[seeda]+indeg[seedb]
    for t in range(3):
        new_y=y.copy()
        for i in range(n):
            # The author's LIM outcome uses one random normal for each node
            # at every time step; draw it even for nodes already adopted.
            e=np.random.normal()
            if y[i]==0:
                exposure=0.
                for j in range(n):
                    if A[i,j]:exposure+=A[i,j]*y[j]
                exposure/=max(1,outdeg[i])
                if alpha+b_degree*seed_degree+b_spill*exposure+e>0:
                    new_y[i]=1
        y=new_y
    return y.sum()/n


@njit
def compute_replicates(A,Din,Dout,Q,norm,ns,alpha,b_degree,beta,N,reps,seed):
    np.random.seed(seed)
    M=len(ns)
    out=np.empty((reps,7))  # truth, DM, HT, Hajek, mean_nom, mean_rand, both_arms
    for r in range(reps):
        order=np.random.permutation(M)
        ys=np.empty(N); ws=np.empty(N,np.uint8)
        wa=np.empty(N);wb=np.empty(N)
        nom_total=0.;rand_total=0.;sum_nom=0.;sum_rand=0.;n_nom=0;n_rand=0
        for j in range(N):
            vi=order[j];n=ns[vi]
            w=np.random.randint(2)
            sa,sb=choose_seeds(Q[vi],n,True)
            ra,rb=choose_seeds(Q[vi],n,False)
            yn=adoption(A[vi],Dout[vi],Din[vi],n,sa,sb,alpha,b_degree,beta)
            yr=adoption(A[vi],Dout[vi],Din[vi],n,ra,rb,alpha,b_degree,beta)
            nom_total+=yn;rand_total+=yr
            ys[j]=yn if w else yr
            ws[j]=w
            a,b=(sa,sb) if w else (ra,rb)
            pa=2*Q[vi,a]*Q[vi,b]/norm[vi]
            pb=2/(n*(n-1))
            den=.5*(pa+pb)
            wa[j]=pa/den;wb[j]=pb/den
            if w:n_nom+=1;sum_nom+=ys[j]
            else:n_rand+=1;sum_rand+=ys[j]
        dm=np.nan
        if n_nom>0 and n_rand>0:dm=sum_nom/n_nom-sum_rand/n_rand
        ht=0.;a_num=0.;b_num=0.;a_den=0.;b_den=0.
        for j in range(N):
            ht+=(wa[j]-wb[j])*ys[j]
            a_num+=wa[j]*ys[j];b_num+=wb[j]*ys[j]
            a_den+=wa[j];b_den+=wb[j]
        ht/=N
        haj=a_num/a_den-b_num/b_den
        out[r,0]=(nom_total-rand_total)/N
        out[r,1]=dm;out[r,2]=ht;out[r,3]=haj
        out[r,4]=nom_total/N;out[r,5]=rand_total/N
        out[r,6]=n_nom
    return out


def main():
    p=argparse.ArgumentParser()
    p.add_argument('--raw',default='data/0422allinforawnet.dta')
    p.add_argument('--replicates',type=int,default=5000)
    p.add_argument('--n-villages',type=int,default=150)
    p.add_argument('--seed',type=int,default=20261102)
    p.add_argument('--degree-effect',type=float,nargs='+',default=list(DEGREE_EFFECTS))
    p.add_argument('--output',type=Path,default=Path('outputs/official_python'))
    args=p.parse_args()
    args.output.mkdir(parents=True,exist_ok=True)
    villages,diag=load_author_villages(args.raw)
    diag.to_csv(args.output/'network_diagnostics.csv',index=False)
    M=len(villages);nmax=max(v[1].shape[0] for v in villages)
    A=np.zeros((M,nmax,nmax),np.uint8)
    Din=np.zeros((M,nmax),np.int64);Dout=np.zeros((M,nmax),np.int64)
    Q=np.zeros((M,nmax));norm=np.zeros(M);ns=np.zeros(M,np.int64)
    for i,(_,a,q,di,do,no) in enumerate(villages):
        n=a.shape[0];ns[i]=n;A[i,:n,:n]=a
        Din[i,:n]=di;Dout[i,:n]=do;Q[i,:n]=q;norm[i]=no
    if args.n_villages>M:raise ValueError('More villages requested than available')
    print(f'Author-filtered villages={M}, mean nodes={ns.mean():.2f}, mean edges={diag.unique_edges.mean():.2f}',flush=True)
    rows=[];summ=[]
    cell=0
    for gamma in args.degree_effect:
        for alpha in ALPHAS:
            for beta in BETAS:
                cell+=1
                arr=compute_replicates(A,Din,Dout,Q,norm,ns,alpha,gamma,beta,args.n_villages,args.replicates,args.seed+cell*101)
                tau=arr[:,0].mean()
                for estimator,index in [('DM',1),('HT',2),('Hajek',3)]:
                    vals=arr[:,index]
                    summ.append(dict(alpha=alpha,beta=beta,b_degree=gamma,estimator=estimator,
                        mean_estimate=np.nanmean(vals),truth=tau,bias=np.nanmean(vals-tau),
                        rmse=np.sqrt(np.nanmean((vals-tau)**2)),
                        onehop=arr[:,4].mean(),random=arr[:,5].mean()))
                    for rep,v in enumerate(vals):
                        rows.append((alpha,beta,gamma,rep,estimator,v,arr[rep,0],tau))
                pd.DataFrame(summ).to_csv(args.output/'summary.csv',index=False)
                print(f'[{cell}/{len(args.degree_effect)*24}] alpha={alpha:g}, beta={beta:g}, b_degree={gamma:g}, truth={tau:.5f}',flush=True)
    pd.DataFrame(rows,columns=['alpha','beta','b_degree','rep','estimator','estimate','replicate_truth','truth']).to_csv(args.output/'replicates.csv',index=False)
    s=pd.DataFrame(summ)
    for gamma in args.degree_effect:
        fig,axes=plt.subplots(2,4,figsize=(16,7),sharex='col')
        for j,alpha in enumerate(ALPHAS):
            x=s[(s.alpha==alpha)&(s.b_degree==gamma)]
            truth=x[x.estimator=='DM'].sort_values('beta')
            axes[0,j].plot(truth.beta,truth.onehop,'o-',label='One-hop')
            axes[0,j].plot(truth.beta,truth.random,'o-',label='Random')
            axes[1,j].plot(truth.beta,truth.truth,'k-',label='Truth')
            for method in ('DM','HT','Hajek'):
                q=x[x.estimator==method].sort_values('beta')
                axes[1,j].plot(q.beta,q.mean_estimate,'o--',label=method,markersize=3)
            axes[0,j].set_title(f'alpha={alpha:g}')
            axes[1,j].set_xlabel('beta')
            for ax in (axes[0,j],axes[1,j]):ax.grid(alpha=.2)
        axes[0,0].set_ylabel('Mean adoption');axes[1,0].set_ylabel('One-hop - Random')
        axes[0,0].legend();axes[1,0].legend()
        fig.suptitle(f'Author-aligned LIM Python replication; b_degree={gamma:g}; N={args.n_villages}')
        fig.tight_layout()
        fig.savefig(args.output/f'figure2_degree_{gamma:g}.png',dpi=160)
        plt.close(fig)
    print('Saved to',args.output.resolve(),flush=True)

if __name__=='__main__':main()
