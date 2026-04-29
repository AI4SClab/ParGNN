
#include "metislib.h"

#define MAXADP 64




void Gvol(ctrl_t *ctrl, graph_t *graph, idx_t niter, 
         real_t ffactor, idx_t omode)
{
  
  idx_t i, ii, iii, j, k, l, pass, nvtxs, nparts, gain; 
  idx_t from, me, to, oldcut, vwgt;
  idx_t *xadj, *adjncy;
  idx_t *where, *pwgts, *perm, *bndptr, *bndind, *minpwgts, *maxpwgts;
  idx_t nmoved, nupd, *vstatus, *updptr, *updind;
  idx_t maxndoms, *safetos=NULL, *nads=NULL, *doms=NULL, **adids=NULL, **adwgts=NULL;
  idx_t *bfslvl=NULL, *bfsind=NULL, *bfsmrk=NULL;
  idx_t bndtype = (omode == OMODE_REFINE ? BNDTYPE_REFINE : BNDTYPE_BALANCE);
  real_t *tpwgts;

  ipq_t *queue;
  idx_t oldvol, xgain;
  idx_t *vmarker, *pmarker, *modind;
  vkrinfo_t *myrinfo;
  vnbr_t *mynbrs;
  WCOREPUSH;


  nvtxs  = graph->nvtxs;
  xadj   = graph->xadj;
  adjncy = graph->adjncy;
  bndptr = graph->bndptr;
  bndind = graph->bndind;
  where  = graph->where;
  pwgts  = graph->pwgts;
  
  nparts = ctrl->nparts;
  tpwgts = ctrl->tpwgts;

 
  minpwgts  = iwspacemalloc(ctrl, nparts);
  maxpwgts  = iwspacemalloc(ctrl, nparts);

  for (i=0; i<nparts; i++) {
    maxpwgts[i]  = ctrl->tpwgts[i]*graph->tvwgt[0]*ctrl->ubfactors[0];
    minpwgts[i]  = ctrl->tpwgts[i]*graph->tvwgt[0]*(1.0/ctrl->ubfactors[0]);
  }

  perm = iwspacemalloc(ctrl, nvtxs);


  
  safetos = iset(nparts, 2, iwspacemalloc(ctrl, nparts));

  if (ctrl->minconn) {
    ComputeSubDomainGraph(ctrl, graph);

    nads    = ctrl->nads;
    adids   = ctrl->adids;
    adwgts  = ctrl->adwgts;
    doms    = iset(nparts, 0, ctrl->pvec1);
  }


  
  vstatus = iset(nvtxs, VPQSTATUS_NOTPRESENT, iwspacemalloc(ctrl, nvtxs));
  updptr  = iset(nvtxs, -1, iwspacemalloc(ctrl, nvtxs));
  updind  = iwspacemalloc(ctrl, nvtxs);

  if (ctrl->contig) {
    
    bfslvl = iset(nvtxs, 0, iwspacemalloc(ctrl, nvtxs));
    bfsind = iwspacemalloc(ctrl, nvtxs);
    bfsmrk = iset(nvtxs, 0, iwspacemalloc(ctrl, nvtxs));
  }

 
  modind  = iwspacemalloc(ctrl, nvtxs);
  vmarker = iset(nvtxs, 0, iwspacemalloc(ctrl, nvtxs));
  pmarker = iset(nparts, -1, iwspacemalloc(ctrl, nparts));

  if (ctrl->dbglvl&METIS_DBG_REFINE) {
     printf("%s: [%6"PRIDX" %6"PRIDX"]-[%6"PRIDX" %6"PRIDX"], Bal: %5.3"PRREAL
         ", Nv-Nb[%6"PRIDX" %6"PRIDX"], Cut: %5"PRIDX", Vol: %5"PRIDX,
         (omode == OMODE_REFINE ? "GRV" : "GBV"),
         pwgts[iargmin(nparts, pwgts,1)], imax(nparts, pwgts,1), minpwgts[0], maxpwgts[0], 
         ComputeLoadImbalance(graph, nparts, ctrl->pijbm), 
         graph->nvtxs, graph->nbnd, graph->mincut, graph->minvol);
     if (ctrl->minconn) 
       printf(", Doms: [%3"PRIDX" %4"PRIDX"]", imax(nparts, nads,1), isum(nparts, nads,1));
     printf("\n");
  }

  queue = ipqCreate(nvtxs);

idx_t nAdp=0;
  
  while (nAdp<=MAXADP) {
    ASSERT(ComputeVolume(graph, where) == graph->minvol);

    if (omode == OMODE_BALANCE) {
      
      for (i=0; i<nparts; i++) {
        if (pwgts[i] > maxpwgts[i])
          break;
      }
      if (i == nparts) 
        break;
    }

    oldcut = graph->mincut;
    oldvol = graph->minvol;
    nupd   = 0;

    if (ctrl->minconn)
      maxndoms = imax(nparts, nads,1);

   
    irandArrayPermute(graph->nbnd, perm, graph->nbnd/4, 1);
    for (ii=0; ii<graph->nbnd; ii++) {
      i = bndind[perm[ii]];
      ipqInsert(queue, i, graph->vkrinfo[i].gv);
      vstatus[i] = VPQSTATUS_PRESENT;
      ListInsert(nupd, updind, updptr, i);
    }

  
    for (nmoved=0, iii=0;;iii++) {
      if ((i = ipqGetTop(queue)) == -1) 
        break;
      vstatus[i] = VPQSTATUS_EXTRACTED;

      myrinfo = graph->vkrinfo+i;
      mynbrs  = ctrl->vnbrpool + myrinfo->inbr;

      from = where[i];
      vwgt = graph->vwgt[i];

     
      if (omode == OMODE_REFINE) {
        if (myrinfo->nid > 0 && pwgts[from]-vwgt < minpwgts[from]) 
          continue;
      }
      else { 
        if (pwgts[from]-vwgt < minpwgts[from]) 
          continue;
      }

      if (ctrl->contig && IsArticulationNode_2(i, xadj, adjncy, where, bfslvl, bfsind, bfsmrk))
        continue;

      if (ctrl->minconn)
        SelectSafeTargetSubdomains(myrinfo, mynbrs, nads, adids, maxndoms, safetos, doms);

      xgain = (myrinfo->nid == 0 && myrinfo->ned > 0 ? graph->vsize[i] : 0);

     
      if (omode == OMODE_REFINE) {
        for (k=myrinfo->nnbrs-1; k>=0; k--) {
          if (!safetos[to=mynbrs[k].pid])
            continue;
           real_t combinedGain = 0;

        
        real_t weightGain = mynbrs[k].gv; 

        
        real_t cutGain = mynbrs[k].ned - myrinfo->nid; 

        
        combinedGain = weightGain + ffactor * cutGain;  
          gain = combinedGain + xgain;
          if (gain >= 0 && pwgts[to]+vwgt <= maxpwgts[to]+ffactor*gain)  
            break;
        }
        if (k < 0)
          continue;  

        for (j=k-1; j>=0; j--) {
          if (!safetos[to=mynbrs[j].pid])
            continue;
          gain = mynbrs[j].gv + xgain;
          if ((mynbrs[j].gv > mynbrs[k].gv && 
               pwgts[to]+vwgt <= maxpwgts[to]+ffactor*gain) 
              ||
              (mynbrs[j].gv == mynbrs[k].gv && 
               mynbrs[j].ned > mynbrs[k].ned &&
               pwgts[to]+vwgt <= maxpwgts[to]) 
              ||
              (mynbrs[j].gv == mynbrs[k].gv && 
               mynbrs[j].ned == mynbrs[k].ned &&
               tpwgts[mynbrs[k].pid]*pwgts[to] < tpwgts[to]*pwgts[mynbrs[k].pid])
             )
            k = j;
        }
        to = mynbrs[k].pid;

        ASSERT(xgain+mynbrs[k].gv >= 0);

        j = 0;
        if (xgain+mynbrs[k].gv > 0 || mynbrs[k].ned-myrinfo->nid > 0)
          j = 1;
        else if (mynbrs[k].ned-myrinfo->nid == 0) {
          if ((iii%2 == 0 && safetos[to] == 2) || 
              pwgts[from] >= maxpwgts[from] || 
              tpwgts[from]*(pwgts[to]+vwgt) < tpwgts[to]*pwgts[from])
            j = 1;
        }
        if (j == 0)
          continue;
      }
      else { 
        for (k=myrinfo->nnbrs-1; k>=0; k--) {
          if (!safetos[to=mynbrs[k].pid])
            continue;
          if (pwgts[to]+vwgt <= maxpwgts[to] || 
              tpwgts[from]*(pwgts[to]+vwgt) <= tpwgts[to]*pwgts[from])  
            break;
        }
        if (k < 0)
          continue;  

        for (j=k-1; j>=0; j--) {
          if (!safetos[to=mynbrs[j].pid])
            continue;
          if (tpwgts[mynbrs[k].pid]*pwgts[to] < tpwgts[to]*pwgts[mynbrs[k].pid])
            k = j;
        }
        to = mynbrs[k].pid;

        if (pwgts[from] < maxpwgts[from] && pwgts[to] > minpwgts[to] && 
            (xgain+mynbrs[k].gv < 0 || 
             (xgain+mynbrs[k].gv == 0 &&  mynbrs[k].ned-myrinfo->nid < 0))
           )
          continue;
      }
          
          
     
      INC_DEC(pwgts[to], pwgts[from], vwgt);
      graph->mincut -= mynbrs[k].ned-myrinfo->nid;
      graph->minvol -= (xgain+mynbrs[k].gv);
      where[i] = to;
      nmoved++;

      IFSET(ctrl->dbglvl, METIS_DBG_MOVEINFO, 
          printf("\t\tMoving %6"PRIDX" from %3"PRIDX" to %3"PRIDX". "
                 "Gain: [%4"PRIDX" %4"PRIDX"]. Cut: %6"PRIDX", Vol: %6"PRIDX"\n", 
              i, from, to, xgain+mynbrs[k].gv, mynbrs[k].ned-myrinfo->nid, 
              graph->mincut, graph->minvol));

     
      if (ctrl->minconn) {
        
        UpdateEdgeSubDomainGraph(ctrl, from, to, myrinfo->nid-mynbrs[k].ned, &maxndoms);

        
        for (j=xadj[i]; j<xadj[i+1]; j++) {
          me = where[adjncy[j]];
          if (me != from && me != to) {
            UpdateEdgeSubDomainGraph(ctrl, from, me, -1, &maxndoms);
            UpdateEdgeSubDomainGraph(ctrl, to, me, 1, &maxndoms);
          }
        }
      }

      
      KWayVolUpdate_2(ctrl, graph, i, from, to, queue, vstatus, &nupd, updptr, 
          updind, bndtype, vmarker, pmarker, modind);

     
    }

nAdp+=nmoved;
 
    for (i=0; i<nupd; i++) {
      ASSERT(updptr[updind[i]] != -1);
      ASSERT(vstatus[updind[i]] != VPQSTATUS_NOTPRESENT);
      vstatus[updind[i]] = VPQSTATUS_NOTPRESENT;
      updptr[updind[i]]  = -1;
    }

    if (ctrl->dbglvl&METIS_DBG_REFINE) {
       printf("\t[%6"PRIDX" %6"PRIDX"], Bal: %5.3"PRREAL", Nb: %6"PRIDX"."
              " Nmoves: %5"PRIDX", Cut: %6"PRIDX", Vol: %6"PRIDX,
              pwgts[iargmin(nparts, pwgts,1)], imax(nparts, pwgts,1),
              ComputeLoadImbalance(graph, nparts, ctrl->pijbm), 
              graph->nbnd, nmoved, graph->mincut, graph->minvol);
       if (ctrl->minconn) 
         printf(", Doms: [%3"PRIDX" %4"PRIDX"]", imax(nparts, nads,1), isum(nparts, nads,1));
       printf("\n");
    }

    if (nmoved == 0 || 
        (omode == OMODE_REFINE && graph->minvol == oldvol && graph->mincut == oldcut))
      break;
  }

  ipqDestroy(queue);

  WCOREPOP;
}



idx_t IsArticulationNode_2(idx_t i, idx_t *xadj, idx_t *adjncy, idx_t *where,
          idx_t *bfslvl, idx_t *bfsind, idx_t *bfsmrk)
{
  idx_t ii, j, k=0, head, tail, nhits, tnhits, from, BFSDEPTH=5;

  from = where[i];

 
  for (tnhits=0, j=xadj[i]; j<xadj[i+1]; j++) {
    if (where[adjncy[j]] == from) {
      ASSERT(bfsmrk[adjncy[j]] == 0);
      ASSERT(bfslvl[adjncy[j]] == 0);
      bfsmrk[k=adjncy[j]] = 1;
      tnhits++;
    }
  }

 
  if (tnhits == 0)
    return 0;
  if (tnhits == 1) {
    bfsmrk[k] = 0;
    return 0;
  }

  ASSERT(bfslvl[i] == 0);
  bfslvl[i] = 1;

  bfsind[0] = k; 
  bfslvl[k] = 1;
  bfsmrk[k] = 0;
  head = 0;
  tail = 1;


  for (nhits=1; head<tail; ) {
    ii = bfsind[head++];
    for (j=xadj[ii]; j<xadj[ii+1]; j++) {
      if (where[k=adjncy[j]] == from) {
        if (bfsmrk[k]) {
          bfsmrk[k] = 0;
          if (++nhits == tnhits)
            break;
        }
        if (bfslvl[k] == 0 && bfslvl[ii] < BFSDEPTH) {
          bfsind[tail++] = k;
          bfslvl[k] = bfslvl[ii]+1;
        }
      }
    }
    if (nhits == tnhits)
      break;
  }

  
  bfslvl[i] = 0;
  for (j=0; j<tail; j++)
    bfslvl[bfsind[j]] = 0;


  
  if (nhits < tnhits) {
    for (j=xadj[i]; j<xadj[i+1]; j++) 
      if (where[adjncy[j]] == from) 
        bfsmrk[adjncy[j]] = 0;
  }

  return (nhits != tnhits);
}



void KWayVolUpdate_2(ctrl_t *ctrl, graph_t *graph, idx_t v, idx_t from, 
         idx_t to, ipq_t *queue, idx_t *vstatus, idx_t *r_nupd, idx_t *updptr, 
         idx_t *updind, idx_t bndtype, idx_t *vmarker, idx_t *pmarker, 
         idx_t *modind)
{
  idx_t i, ii, iii, j, jj, k, kk, l, u, nmod, other, me, myidx; 
  idx_t *xadj, *vsize, *adjncy, *where;
  vkrinfo_t *myrinfo, *orinfo;
  vnbr_t *mynbrs, *onbrs;

  xadj   = graph->xadj;
  adjncy = graph->adjncy;
  vsize  = graph->vsize;
  where  = graph->where;

  myrinfo = graph->vkrinfo+v;
  mynbrs  = ctrl->vnbrpool + myrinfo->inbr;


  
  for (k=0; k<myrinfo->nnbrs; k++)
    pmarker[mynbrs[k].pid] = k;
  pmarker[from] = k;

  myidx = pmarker[to];  

  for (j=xadj[v]; j<xadj[v+1]; j++) {
    ii     = adjncy[j];
    other  = where[ii];
    orinfo = graph->vkrinfo+ii;
    onbrs  = ctrl->vnbrpool + orinfo->inbr;

    if (other == from) {
      for (k=0; k<orinfo->nnbrs; k++) {
        if (pmarker[onbrs[k].pid] == -1) 
          onbrs[k].gv += vsize[v];
      }
    }
    else {
      ASSERT(pmarker[other] != -1);

      if (mynbrs[pmarker[other]].ned > 1) {
        for (k=0; k<orinfo->nnbrs; k++) {
          if (pmarker[onbrs[k].pid] == -1) 
            onbrs[k].gv += vsize[v];
        }
      }
      else { 
        for (k=0; k<orinfo->nnbrs; k++) {
          if (pmarker[onbrs[k].pid] != -1) 
            onbrs[k].gv -= vsize[v];
        }
      }
    }
  }

  for (k=0; k<myrinfo->nnbrs; k++)
    pmarker[mynbrs[k].pid] = -1;
  pmarker[from] = -1;


  
  if (myidx == -1) {
    myidx = myrinfo->nnbrs++;
    ASSERT(myidx < xadj[v+1]-xadj[v]);
    mynbrs[myidx].ned = 0;
  }
  myrinfo->ned += myrinfo->nid-mynbrs[myidx].ned;
  SWAP(myrinfo->nid, mynbrs[myidx].ned, j);
  if (mynbrs[myidx].ned == 0) 
    mynbrs[myidx] = mynbrs[--myrinfo->nnbrs];
  else
    mynbrs[myidx].pid = from;


  
  vmarker[v] = 1;
  modind[0]  = v;
  nmod       = 1;
  for (j=xadj[v]; j<xadj[v+1]; j++) {
    ii = adjncy[j];
    me = where[ii];

    if (!vmarker[ii]) {  
      vmarker[ii] = 2;
      modind[nmod++] = ii;
    }

    myrinfo = graph->vkrinfo+ii;
    if (myrinfo->inbr == -1) 
      myrinfo->inbr = vnbrpoolGetNext(ctrl, xadj[ii+1]-xadj[ii]);
    mynbrs = ctrl->vnbrpool + myrinfo->inbr;

    if (me == from) {
      INC_DEC(myrinfo->ned, myrinfo->nid, 1);
    } 
    else if (me == to) {
      INC_DEC(myrinfo->nid, myrinfo->ned, 1);
    }

   
    if (me != from) {
      for (k=0; k<myrinfo->nnbrs; k++) {
        if (mynbrs[k].pid == from) {
          if (mynbrs[k].ned == 1) {
            mynbrs[k] = mynbrs[--myrinfo->nnbrs];
            vmarker[ii] = 1;  

           
            for (jj=xadj[ii]; jj<xadj[ii+1]; jj++) {
              u      = adjncy[jj];
              other  = where[u];
              orinfo = graph->vkrinfo+u;
              onbrs  = ctrl->vnbrpool + orinfo->inbr;

              for (kk=0; kk<orinfo->nnbrs; kk++) {
                if (onbrs[kk].pid == from) {
                  onbrs[kk].gv -= vsize[ii];
                  if (!vmarker[u]) { 
                    vmarker[u]      = 2;
                    modind[nmod++] = u;
                  }
                  break;
                }
              }
            }
          }
          else {
            mynbrs[k].ned--;

            
            if (mynbrs[k].ned == 1) {
             
              for (jj=xadj[ii]; jj<xadj[ii+1]; jj++) {
                u     = adjncy[jj];
                other = where[u];

                if (other == from) {
                  orinfo = graph->vkrinfo+u;
                  onbrs  = ctrl->vnbrpool + orinfo->inbr;

                 
                  for (kk=0; kk<orinfo->nnbrs; kk++) 
                    onbrs[kk].gv += vsize[ii];

                  if (!vmarker[u]) { 
                    vmarker[u]     = 2;
                    modind[nmod++] = u;
                  }
                  break;  
                }
              }
            }
          }
          break; 
        }
      }
    }


    
    if (me != to) {
      for (k=0; k<myrinfo->nnbrs; k++) {
        if (mynbrs[k].pid == to) {
          mynbrs[k].ned++;

         
          if (mynbrs[k].ned == 2) {
        
            for (jj=xadj[ii]; jj<xadj[ii+1]; jj++) {
              u     = adjncy[jj];
              other = where[u];

              if (u != v && other == to) {
                orinfo = graph->vkrinfo+u;
                onbrs  = ctrl->vnbrpool + orinfo->inbr;
                for (kk=0; kk<orinfo->nnbrs; kk++) 
                  onbrs[kk].gv -= vsize[ii];

                if (!vmarker[u]) { 
                  vmarker[u]      = 2;
                  modind[nmod++] = u;
                }
                break;  
              }
            }
          }
          break;
        }
      }

      if (k == myrinfo->nnbrs) {
        mynbrs[myrinfo->nnbrs].pid   = to;
        mynbrs[myrinfo->nnbrs++].ned = 1;
        vmarker[ii] = 1;  
        for (jj=xadj[ii]; jj<xadj[ii+1]; jj++) {
          u      = adjncy[jj];
          other  = where[u];
          orinfo = graph->vkrinfo+u;
          onbrs  = ctrl->vnbrpool + orinfo->inbr;

          for (kk=0; kk<orinfo->nnbrs; kk++) {
            if (onbrs[kk].pid == to) {
              onbrs[kk].gv += vsize[ii];
              if (!vmarker[u]) { /* Need to update boundary etc */
                vmarker[u] = 2;
                modind[nmod++] = u;
              }
              break;
            }
          }
        }
      }
    }

    ASSERT(myrinfo->nnbrs <= xadj[ii+1]-xadj[ii]);
  }


  
  myrinfo = graph->vkrinfo+v;
  mynbrs  = ctrl->vnbrpool + myrinfo->inbr;
  for (k=0; k<myrinfo->nnbrs; k++)
    pmarker[mynbrs[k].pid] = k;
  pmarker[to] = k;

  for (j=xadj[v]; j<xadj[v+1]; j++) {
    ii     = adjncy[j];
    other  = where[ii];
    orinfo = graph->vkrinfo+ii;
    onbrs  = ctrl->vnbrpool + orinfo->inbr;

    if (other == to) {
      for (k=0; k<orinfo->nnbrs; k++) {
        if (pmarker[onbrs[k].pid] == -1) 
          onbrs[k].gv -= vsize[v];
      }
    }
    else {
      ASSERT(pmarker[other] != -1);

      if (mynbrs[pmarker[other]].ned > 1) {
        for (k=0; k<orinfo->nnbrs; k++) {
          if (pmarker[onbrs[k].pid] == -1) 
            onbrs[k].gv -= vsize[v];
        }
      }
      else { 
        for (k=0; k<orinfo->nnbrs; k++) {
          if (pmarker[onbrs[k].pid] != -1) 
            onbrs[k].gv += vsize[v];
        }
      }
    }
  }
  for (k=0; k<myrinfo->nnbrs; k++)
    pmarker[mynbrs[k].pid] = -1;
  pmarker[to] = -1;


 
  for (iii=0; iii<nmod; iii++) {
    i  = modind[iii];
    me = where[i];

    myrinfo = graph->vkrinfo+i;
    mynbrs  = ctrl->vnbrpool + myrinfo->inbr;

    if (vmarker[i] == 1) {  
      for (k=0; k<myrinfo->nnbrs; k++) 
        mynbrs[k].gv = 0;

      for (j=xadj[i]; j<xadj[i+1]; j++) {
        ii     = adjncy[j];
        other  = where[ii];
        orinfo = graph->vkrinfo+ii;
        onbrs  = ctrl->vnbrpool + orinfo->inbr;

        for (kk=0; kk<orinfo->nnbrs; kk++) 
          pmarker[onbrs[kk].pid] = kk;
        pmarker[other] = 1;

        if (me == other) {
          
          for (k=0; k<myrinfo->nnbrs; k++) {
            if (pmarker[mynbrs[k].pid] == -1)
              mynbrs[k].gv -= vsize[ii];
          }
        }
        else {
          ASSERT(pmarker[me] != -1);

        
          if (onbrs[pmarker[me]].ned == 1) { 
        
            for (k=0; k<myrinfo->nnbrs; k++) {
              if (pmarker[mynbrs[k].pid] != -1) 
                mynbrs[k].gv += vsize[ii];
            }
          }
          else {
            
            for (k=0; k<myrinfo->nnbrs; k++) {
              if (pmarker[mynbrs[k].pid] == -1) 
                mynbrs[k].gv -= vsize[ii];
            }
          }
        }

        for (kk=0; kk<orinfo->nnbrs; kk++) 
          pmarker[onbrs[kk].pid] = -1;
        pmarker[other] = -1;
  
      }
    }

   
    myrinfo->gv = IDX_MIN;
    for (k=0; k<myrinfo->nnbrs; k++) {
      if (mynbrs[k].gv > myrinfo->gv)
        myrinfo->gv = mynbrs[k].gv;
    }

   
    if (myrinfo->ned > 0 && myrinfo->nid == 0)
      myrinfo->gv += vsize[i];


    
    if (bndtype == BNDTYPE_REFINE) {
      if (myrinfo->gv >= 0 && graph->bndptr[i] == -1)
        BNDInsert(graph->nbnd, graph->bndind, graph->bndptr, i);

      if (myrinfo->gv < 0 && graph->bndptr[i] != -1)
        BNDDelete(graph->nbnd, graph->bndind, graph->bndptr, i);
    }
    else {
      if (myrinfo->ned > 0 && graph->bndptr[i] == -1)
        BNDInsert(graph->nbnd, graph->bndind, graph->bndptr, i);

      if (myrinfo->ned == 0 && graph->bndptr[i] != -1)
        BNDDelete(graph->nbnd, graph->bndind, graph->bndptr, i);
    }


   
    if (queue != NULL) {
      if (vstatus[i] != VPQSTATUS_EXTRACTED) {
        if (graph->bndptr[i] != -1) { 
          if (vstatus[i] == VPQSTATUS_PRESENT) {
            ipqUpdate(queue, i, myrinfo->gv);
          }
          else {
            ipqInsert(queue, i, myrinfo->gv);
            vstatus[i] = VPQSTATUS_PRESENT;
            ListInsert(*r_nupd, updind, updptr, i);
          }
        }
        else { 
          if (vstatus[i] == VPQSTATUS_PRESENT) {
            ipqDelete(queue, i);
            vstatus[i] = VPQSTATUS_NOTPRESENT;
            ListDelete(*r_nupd, updind, updptr, i);
          }
        }
      }
    }
  
    vmarker[i] = 0;
  }
}





void xinGreedy_KWayCutOptimize(ctrl_t *ctrl, graph_t *graph, idx_t niter, 
         real_t ffactor, idx_t omode)
{
 
  idx_t i, ii, iii, j, k, l, pass, nvtxs, nparts, gain; 
  idx_t from, me, to, oldcut, vwgt;
  idx_t *xadj, *adjncy, *adjwgt;
  idx_t *where, *pwgts, *perm, *bndptr, *bndind, *minpwgts, *maxpwgts;
  idx_t nmoved, nupd, *vstatus, *updptr, *updind;
  idx_t maxndoms, *safetos=NULL, *nads=NULL, *doms=NULL, **adids=NULL, **adwgts=NULL;
  idx_t *bfslvl=NULL, *bfsind=NULL, *bfsmrk=NULL;
  idx_t bndtype = (omode == OMODE_REFINE ? BNDTYPE_REFINE : BNDTYPE_BALANCE);
  real_t *tpwgts, ubfactor;


  idx_t nbnd, oldnnbrs;
  rpq_t *queue;
  real_t rgain;
  ckrinfo_t *myrinfo;
  cnbr_t *mynbrs;

  ffactor = 0.0;
  WCOREPUSH;
    //*********************************************** */  
    // 在函数开始时,为了后面比较重划分后和原始划分的区别

    idx_t *initial_where = imalloc(graph->nvtxs, "initial_where");
    icopy(graph->nvtxs, graph->where, initial_where);
    //*********************************************** */ 
  nvtxs  = graph->nvtxs;
  xadj   = graph->xadj;
  adjncy = graph->adjncy;
  adjwgt = graph->adjwgt;

  bndind = graph->bndind;
  bndptr = graph->bndptr;

  where = graph->where;
  pwgts = graph->pwgts;
  
  nparts = ctrl->nparts;
  tpwgts = ctrl->tpwgts;

  
  minpwgts = iwspacemalloc(ctrl, nparts);
  maxpwgts = iwspacemalloc(ctrl, nparts);

  if (omode == OMODE_BALANCE)
    ubfactor = ctrl->ubfactors[0];
  else
    ubfactor = gk_max(ctrl->ubfactors[0], ComputeLoadImbalance(graph, nparts, ctrl->pijbm));

  for (i=0; i<nparts; i++) {
    maxpwgts[i] = tpwgts[i]*graph->tvwgt[0]*ubfactor;
    minpwgts[i] = tpwgts[i]*graph->tvwgt[0]*(1.0/ubfactor);
  }

  perm = iwspacemalloc(ctrl, nvtxs);


  
  safetos = iset(nparts, 2, iwspacemalloc(ctrl, nparts));

  if (ctrl->minconn) {
    ComputeSubDomainGraph(ctrl, graph);

    nads    = ctrl->nads;
    adids   = ctrl->adids;
    adwgts  = ctrl->adwgts;
    doms    = iset(nparts, 0, ctrl->pvec1);
  }


  
  vstatus = iset(nvtxs, VPQSTATUS_NOTPRESENT, iwspacemalloc(ctrl, nvtxs));
  updptr  = iset(nvtxs, -1, iwspacemalloc(ctrl, nvtxs));
  updind  = iwspacemalloc(ctrl, nvtxs);

  if (ctrl->contig) {
    
    bfslvl = iset(nvtxs, 0, iwspacemalloc(ctrl, nvtxs));
    bfsind = iwspacemalloc(ctrl, nvtxs);
    bfsmrk = iset(nvtxs, 0, iwspacemalloc(ctrl, nvtxs));
  }

  if (ctrl->dbglvl&METIS_DBG_REFINE) {
     printf("%s: [%6"PRIDX" %6"PRIDX"]-[%6"PRIDX" %6"PRIDX"], Bal: %5.3"PRREAL"," 
            " Nv-Nb[%6"PRIDX" %6"PRIDX"], Cut: %6"PRIDX,
            (omode == OMODE_REFINE ? "GRC" : "GBC"),
            pwgts[iargmin(nparts, pwgts,1)], imax(nparts, pwgts,1), minpwgts[0], maxpwgts[0], 
            ComputeLoadImbalance(graph, nparts, ctrl->pijbm), 
            graph->nvtxs, graph->nbnd, graph->mincut);
     if (ctrl->minconn) 
       printf(", Doms: [%3"PRIDX" %4"PRIDX"]", imax(nparts, nads,1), isum(nparts, nads,1));
     printf("\n");
  }

  queue = rpqCreate(nvtxs);

 
   idx_t nAdp = 0;
   ctrl->dbglvl |= METIS_DBG_REFINE;
  while (nAdp < MAXADP) {
    
    ASSERT(ComputeCut(graph, where) == graph->mincut);
    if (omode == OMODE_REFINE)
      ASSERT(CheckBnd2(graph));

    if (omode == OMODE_BALANCE) {
     
      for (i=0; i<nparts; i++) {
        if (pwgts[i] > maxpwgts[i] || pwgts[i] < minpwgts[i])
          break;
      }
      if (i == nparts) 
        break;
    }

    oldcut = graph->mincut;
    nbnd   = graph->nbnd;
    nupd   = 0;

    if (ctrl->minconn)
      maxndoms = imax(nparts, nads,1);

    
    irandArrayPermute(nbnd, perm, nbnd/4, 1);
    for (ii=0; ii<nbnd; ii++) {
      i = bndind[perm[ii]];
      rgain = (graph->ckrinfo[i].nnbrs > 0 ? 
               1.0*graph->ckrinfo[i].ed/sqrt(graph->ckrinfo[i].nnbrs) : 0.0) 
               - graph->ckrinfo[i].id;
      rpqInsert(queue, i, rgain);
      vstatus[i] = VPQSTATUS_PRESENT;
      ListInsert(nupd, updind, updptr, i);
    }

    
    for (nmoved=0, iii=0;;iii++) {
      if ((i = rpqGetTop(queue)) == -1) 
        break;
      vstatus[i] = VPQSTATUS_EXTRACTED;

      myrinfo = graph->ckrinfo+i;
      mynbrs  = ctrl->cnbrpool + myrinfo->inbr;

      from = where[i];
      vwgt = graph->vwgt[i];

#ifdef XXX
      
      if (omode == OMODE_REFINE) {
        if (myrinfo->id > 0 && pwgts[from]-vwgt < minpwgts[from]) 
          continue;   
      }
      else { 
        if (pwgts[from]-vwgt < minpwgts[from]) 
          continue;   
      }
#endif

      if (ctrl->contig && IsArticulationNode_2(i, xadj, adjncy, where, bfslvl, bfsind, bfsmrk))
        continue;

      if (ctrl->minconn)
        SelectSafeTargetSubdomains(myrinfo, mynbrs, nads, adids, maxndoms, safetos, doms);

     
      if (omode == OMODE_REFINE) {
        for (k=myrinfo->nnbrs-1; k>=0; k--) {
          if (!safetos[to=mynbrs[k].pid])
            continue;
          if (((mynbrs[k].ed > myrinfo->id) && 
               ((pwgts[from]-vwgt >= minpwgts[from]) ||
                (tpwgts[from]*pwgts[to] < tpwgts[to]*(pwgts[from]-vwgt))) &&
               ((pwgts[to]+vwgt <= maxpwgts[to]) || 
                (tpwgts[from]*pwgts[to] < tpwgts[to]*(pwgts[from]-vwgt)))
              ) ||
              ((mynbrs[k].ed == myrinfo->id) && 
               (tpwgts[from]*pwgts[to] < tpwgts[to]*(pwgts[from]-vwgt)))
             )
            break;
        }
        if (k < 0)
          continue;  

        for (j=k-1; j>=0; j--) {
          if (!safetos[to=mynbrs[j].pid])
            continue;
          if (((mynbrs[j].ed > mynbrs[k].ed) && 
              ((pwgts[from]-vwgt >= minpwgts[from]) ||
               (tpwgts[from]*pwgts[to] < tpwgts[to]*(pwgts[from]-vwgt))) &&
              ((pwgts[to]+vwgt <= maxpwgts[to]) ||
               (tpwgts[from]*pwgts[to] < tpwgts[to]*(pwgts[from]-vwgt))) 
              ) ||
              ((mynbrs[j].ed == mynbrs[k].ed) && 
               (tpwgts[mynbrs[k].pid]*pwgts[to] < tpwgts[to]*pwgts[mynbrs[k].pid]))
             )
            k = j;
        }

        to = mynbrs[k].pid;

        gain = mynbrs[k].ed-myrinfo->id;
       
      }
      else {  
        for (k=myrinfo->nnbrs-1; k>=0; k--) {
          if (!safetos[to=mynbrs[k].pid])
            continue;
          
          if (from >= nparts || tpwgts[from]*pwgts[to] < tpwgts[to]*(pwgts[from]-vwgt))
            break;
        }
        if (k < 0)
          continue;  

        for (j=k-1; j>=0; j--) {
          if (!safetos[to=mynbrs[j].pid])
            continue;
          if (tpwgts[mynbrs[k].pid]*pwgts[to] < tpwgts[to]*pwgts[mynbrs[k].pid]) 
            k = j;
        }

        to = mynbrs[k].pid;

       
      }


      
      graph->mincut -= mynbrs[k].ed-myrinfo->id;
      nmoved++;

      IFSET(ctrl->dbglvl, METIS_DBG_MOVEINFO, 
          printf("\t\tMoving %6"PRIDX" from %3"PRIDX"/%"PRIDX" to %3"PRIDX"/%"PRIDX" [%6"PRIDX" %6"PRIDX"]. Gain: %4"PRIDX". Cut: %6"PRIDX"\n", 
              i, from, safetos[from], to, safetos[to], pwgts[from], pwgts[to], mynbrs[k].ed-myrinfo->id, graph->mincut));

      
      if (ctrl->minconn) {
        
        UpdateEdgeSubDomainGraph(ctrl, from, to, myrinfo->id-mynbrs[k].ed, &maxndoms);

      
        for (j=xadj[i]; j<xadj[i+1]; j++) {
          me = where[adjncy[j]];
          if (me != from && me != to) {
            UpdateEdgeSubDomainGraph(ctrl, from, me, -adjwgt[j], &maxndoms);
            UpdateEdgeSubDomainGraph(ctrl, to, me, adjwgt[j], &maxndoms);
          }
        }
      }

    
      INC_DEC(pwgts[to], pwgts[from], vwgt);
      UpdateMovedVertexInfoAndBND(i, from, k, to, myrinfo, mynbrs, where, nbnd, 
          bndptr, bndind, bndtype);
      
      
      for (j=xadj[i]; j<xadj[i+1]; j++) {
        ii = adjncy[j];
        me = where[ii];
        myrinfo = graph->ckrinfo+ii;

        oldnnbrs = myrinfo->nnbrs;

        UpdateAdjacentVertexInfoAndBND(ctrl, ii, xadj[ii+1]-xadj[ii], me, 
            from, to, myrinfo, adjwgt[j], nbnd, bndptr, bndind, bndtype);

        UpdateQueueInfo(queue, vstatus, ii, me, from, to, myrinfo, oldnnbrs, 
            nupd, updptr, updind, bndtype);

        ASSERT(myrinfo->nnbrs <= xadj[ii+1]-xadj[ii]);
      }

    }
     // 迭代结束后，更新微调点数的计数器
    nAdp += nmoved;

    graph->nbnd = nbnd;

   
    for (i=0; i<nupd; i++) {
      ASSERT(updptr[updind[i]] != -1);
      ASSERT(vstatus[updind[i]] != VPQSTATUS_NOTPRESENT);
      vstatus[updind[i]] = VPQSTATUS_NOTPRESENT;
      updptr[updind[i]]  = -1;
    }

    

    
 
 

    // 记得释放 initial_where 内存
    gk_free((void **)&initial_where, LTERM);
    break; 
    }
  rpqDestroy(queue);

  WCOREPOP;
}
   


// }



void GKO(ctrl_t *ctrl, graph_t *graph, idx_t niter, 
         real_t ffactor, idx_t omode)
{
  switch (ctrl->objtype) {
    case METIS_OBJTYPE_CUT:
      if (graph->ncon == 1)
        {xinGreedy_KWayCutOptimize(ctrl, graph, niter, ffactor, omode);}
      

    case METIS_OBJTYPE_VOL:
      if (graph->ncon == 1)
       { Gvol(ctrl, graph, niter, ffactor, omode);}
      else
         {Gvol(ctrl, graph, niter, ffactor, omode);}
      break;

    default:
      gk_errexit(SIGERR, "Unknown objtype of %d\n", ctrl->objtype);
  }
}




void AllocateKWayPartitionMemory_2(ctrl_t *ctrl, graph_t *graph)
{

  graph->pwgts  = imalloc(ctrl->nparts*graph->ncon, "AllocateKWayPartitionMemory_2: pwgts");
  graph->where  = imalloc(graph->nvtxs,  "AllocateKWayPartitionMemory_2: where");
  graph->bndptr = imalloc(graph->nvtxs,  "AllocateKWayPartitionMemory_2: bndptr");
  graph->bndind = imalloc(graph->nvtxs,  "AllocateKWayPartitionMemory_2: bndind");

  switch (ctrl->objtype) {
    case METIS_OBJTYPE_CUT:
      graph->ckrinfo  = (ckrinfo_t *)gk_malloc(graph->nvtxs*sizeof(ckrinfo_t), 
                          "AllocateKWayPartitionMemory_2: ckrinfo");
      break;

    case METIS_OBJTYPE_VOL:
      graph->vkrinfo = (vkrinfo_t *)gk_malloc(graph->nvtxs*sizeof(vkrinfo_t), 
                          "AllocateKWayVolPartitionMemory: vkrinfo");

    
      graph->ckrinfo = (ckrinfo_t *)graph->vkrinfo;
      break;

    default:
      gk_errexit(SIGERR, "Unknown objtype of %d\n", ctrl->objtype);
  }

}



void ComputeKWayPartitionParams_2(ctrl_t *ctrl, graph_t *graph)
{
  idx_t i, j, k, l, nvtxs, ncon, nparts, nbnd, mincut, me, other;
  idx_t *xadj, *vwgt, *adjncy, *adjwgt, *pwgts, *where, *bndind, *bndptr;

  nparts = ctrl->nparts;

  nvtxs  = graph->nvtxs;
  ncon   = graph->ncon;
  xadj   = graph->xadj;
  vwgt   = graph->vwgt;
  adjncy = graph->adjncy;
  adjwgt = graph->adjwgt;

  where  = graph->where;
  pwgts  = iset(nparts*ncon, 0, graph->pwgts);
  bndind = graph->bndind;
  bndptr = iset(nvtxs, -1, graph->bndptr);

  nbnd = mincut = 0;

 
  if (ncon == 1) {
    for (i=0; i<nvtxs; i++) {
      ASSERT(where[i] >= 0 && where[i] < nparts);
      pwgts[where[i]] += vwgt[i];
    }
  }
  else {
    for (i=0; i<nvtxs; i++) {
      me = where[i];
      for (j=0; j<ncon; j++)
        pwgts[me*ncon+j] += vwgt[i*ncon+j];
    }
  }


  switch (ctrl->objtype) {
    case METIS_OBJTYPE_CUT:
      {
        ckrinfo_t *myrinfo;
        cnbr_t *mynbrs;

        memset(graph->ckrinfo, 0, sizeof(ckrinfo_t)*nvtxs);
        cnbrpoolReset(ctrl);

        for (i=0; i<nvtxs; i++) {
          me      = where[i];
          myrinfo = graph->ckrinfo+i;

          for (j=xadj[i]; j<xadj[i+1]; j++) {
            if (me == where[adjncy[j]])
              myrinfo->id += adjwgt[j];
            else
              myrinfo->ed += adjwgt[j];
          }

          
          if (myrinfo->ed > 0) {
            mincut += myrinfo->ed;

            myrinfo->inbr = cnbrpoolGetNext(ctrl, xadj[i+1]-xadj[i]);
            mynbrs        = ctrl->cnbrpool + myrinfo->inbr;

            for (j=xadj[i]; j<xadj[i+1]; j++) {
              other = where[adjncy[j]];
              if (me != other) {
                for (k=0; k<myrinfo->nnbrs; k++) {
                  if (mynbrs[k].pid == other) {
                    mynbrs[k].ed += adjwgt[j];
                    break;
                  }
                }
                if (k == myrinfo->nnbrs) {
                  mynbrs[k].pid = other;
                  mynbrs[k].ed  = adjwgt[j];
                  myrinfo->nnbrs++;
                }
              }
            }

            ASSERT(myrinfo->nnbrs <= xadj[i+1]-xadj[i]);

            
            if (myrinfo->ed-myrinfo->id >= 0)
              BNDInsert(nbnd, bndind, bndptr, i);
          }
          else {
            myrinfo->inbr = -1;
          }
        }

        graph->mincut = mincut/2;
        graph->nbnd   = nbnd;

      }
      ASSERT(CheckBnd2(graph));
      break;

    case METIS_OBJTYPE_VOL:
      {
        vkrinfo_t *myrinfo;
        vnbr_t *mynbrs;

        memset(graph->vkrinfo, 0, sizeof(vkrinfo_t)*nvtxs);
        vnbrpoolReset(ctrl);

        /* Compute now the id/ed degrees */
        for (i=0; i<nvtxs; i++) {
          me      = where[i];
          myrinfo = graph->vkrinfo+i;
      
          for (j=xadj[i]; j<xadj[i+1]; j++) {
            if (me == where[adjncy[j]]) 
              myrinfo->nid++;
            else 
              myrinfo->ned++;
          }
      
      
          if (myrinfo->ned > 0) { 
            mincut += myrinfo->ned;

            myrinfo->inbr = vnbrpoolGetNext(ctrl, xadj[i+1]-xadj[i]);
            mynbrs        = ctrl->vnbrpool + myrinfo->inbr;

            for (j=xadj[i]; j<xadj[i+1]; j++) {
              other = where[adjncy[j]];
              if (me != other) {
                for (k=0; k<myrinfo->nnbrs; k++) {
                  if (mynbrs[k].pid == other) {
                    mynbrs[k].ned++;
                    break;
                  }
                }
                if (k == myrinfo->nnbrs) {
                  mynbrs[k].gv  = 0;
                  mynbrs[k].pid = other;
                  mynbrs[k].ned = 1;
                  myrinfo->nnbrs++;
                }
              }
            }
            ASSERT(myrinfo->nnbrs <= xadj[i+1]-xadj[i]);
          }
          else {
            myrinfo->inbr = -1;
          }
        }
        graph->mincut = mincut/2;
      
        ComputeKWayVolGains_2(ctrl, graph);
      }
      ASSERT(graph->minvol == ComputeVolume(graph, graph->where));
      break;
    default:
      gk_errexit(SIGERR, "Unknown objtype of %d\n", ctrl->objtype);
  }

}



void ProjectKWayPartition_2(ctrl_t *ctrl, graph_t *graph)
{
  idx_t i, j, k, nvtxs, nbnd, nparts, me, other, istart, iend, tid, ted;
  idx_t *xadj, *adjncy, *adjwgt;
  idx_t *cmap, *where, *bndptr, *bndind, *cwhere, *htable;
  graph_t *cgraph;
  int dropedges;

  WCOREPUSH;

  dropedges = ctrl->dropedges;

  nparts = ctrl->nparts;

  cgraph = graph->coarser;
  cwhere = cgraph->where;

  if (ctrl->objtype == METIS_OBJTYPE_CUT) {
    ASSERT(CheckBnd2(cgraph));
  }
  else {
    ASSERT(cgraph->minvol == ComputeVolume(cgraph, cgraph->where));
  }


  FreeSData(cgraph);

  nvtxs   = graph->nvtxs;
  cmap    = graph->cmap;
  xadj    = graph->xadj;
  adjncy  = graph->adjncy;
  adjwgt  = graph->adjwgt;

  AllocateKWayPartitionMemory_2(ctrl, graph);

  where  = graph->where;
  bndind = graph->bndind;
  bndptr = iset(nvtxs, -1, graph->bndptr);

  htable = iset(nparts, -1, iwspacemalloc(ctrl, nparts));

  
  switch (ctrl->objtype) {
    case METIS_OBJTYPE_CUT:
      {
        ckrinfo_t *myrinfo;
        cnbr_t *mynbrs;

        
        for (i=0; i<nvtxs; i++) {
          k        = cmap[i];
          where[i] = cwhere[k];
          cmap[i]  = (dropedges ? 1 : cgraph->ckrinfo[k].ed);  
        }

        memset(graph->ckrinfo, 0, sizeof(ckrinfo_t)*nvtxs);
        cnbrpoolReset(ctrl);

        for (nbnd=0, i=0; i<nvtxs; i++) {
          istart = xadj[i];
          iend   = xadj[i+1];

          myrinfo = graph->ckrinfo+i;

          if (cmap[i] == 0) { 
            for (tid=0, j=istart; j<iend; j++) 
              tid += adjwgt[j];

            myrinfo->id   = tid;
            myrinfo->inbr = -1;
          }
          else { 
            myrinfo->inbr = cnbrpoolGetNext(ctrl, iend-istart);
            mynbrs        = ctrl->cnbrpool + myrinfo->inbr;

            me = where[i];
            for (tid=0, ted=0, j=istart; j<iend; j++) {
              other = where[adjncy[j]];
              if (me == other) {
                tid += adjwgt[j];
              }
              else {
                ted += adjwgt[j];
                if ((k = htable[other]) == -1) {
                  htable[other]               = myrinfo->nnbrs;
                  mynbrs[myrinfo->nnbrs].pid  = other;
                  mynbrs[myrinfo->nnbrs++].ed = adjwgt[j];
                }
                else {
                  mynbrs[k].ed += adjwgt[j];
                }
              }
            }
            myrinfo->id = tid;
            myrinfo->ed = ted;
      
            
            if (ted == 0) { 
              ctrl->nbrpoolcpos -= gk_min(nparts, iend-istart);
              myrinfo->inbr      = -1;
            }
            else {
              if (ted-tid >= 0) 
                BNDInsert(nbnd, bndind, bndptr, i); 
      
              for (j=0; j<myrinfo->nnbrs; j++)
                htable[mynbrs[j].pid] = -1;
            }
          }
        }
      
        graph->nbnd = nbnd;
      }
      ASSERT(CheckBnd2(graph));
      break;

    case METIS_OBJTYPE_VOL:
      {
        vkrinfo_t *myrinfo;
        vnbr_t *mynbrs;

       
        for (i=0; i<nvtxs; i++) {
          k        = cmap[i];
          where[i] = cwhere[k];
          cmap[i]  = (dropedges ? 1 : cgraph->vkrinfo[k].ned);  
        }

        memset(graph->vkrinfo, 0, sizeof(vkrinfo_t)*nvtxs);
        vnbrpoolReset(ctrl);

        for (i=0; i<nvtxs; i++) {
          istart = xadj[i];
          iend   = xadj[i+1];
          myrinfo = graph->vkrinfo+i;

          if (cmap[i] == 0) { 
            myrinfo->nid  = iend-istart;
            myrinfo->inbr = -1;
          }
          else { 
            myrinfo->inbr = vnbrpoolGetNext(ctrl, iend-istart);
            mynbrs        = ctrl->vnbrpool + myrinfo->inbr;

            me = where[i];
            for (tid=0, ted=0, j=istart; j<iend; j++) {
              other = where[adjncy[j]];
              if (me == other) {
                tid++;
              }
              else {
                ted++;
                if ((k = htable[other]) == -1) {
                  htable[other]                = myrinfo->nnbrs;
                  mynbrs[myrinfo->nnbrs].gv    = 0;
                  mynbrs[myrinfo->nnbrs].pid   = other;
                  mynbrs[myrinfo->nnbrs++].ned = 1;
                }
                else {
                  mynbrs[k].ned++;
                }
              }
            }
            myrinfo->nid = tid;
            myrinfo->ned = ted;
      
           
            if (ted == 0) { 
              ctrl->nbrpoolcpos -= gk_min(nparts, iend-istart);
              myrinfo->inbr = -1;
            }
            else {
              for (j=0; j<myrinfo->nnbrs; j++)
                htable[mynbrs[j].pid] = -1;
            }
          }
        }
      
        ComputeKWayVolGains_2(ctrl, graph);

        ASSERT(graph->minvol == ComputeVolume(graph, graph->where));
      }
      break;

    default:
      gk_errexit(SIGERR, "Unknown objtype of %d\n", ctrl->objtype);
  }

  graph->mincut = (dropedges ? ComputeCut(graph, where) : cgraph->mincut);
  icopy(nparts*graph->ncon, cgraph->pwgts, graph->pwgts);

  FreeGraph(&graph->coarser);

  WCOREPOP;
}



void ComputeKWayBoundary_2(ctrl_t *ctrl, graph_t *graph, idx_t bndtype)
{
  idx_t i, nvtxs, nbnd;
  idx_t *bndind, *bndptr;

  nvtxs  = graph->nvtxs;
  bndind = graph->bndind;
  bndptr = iset(nvtxs, -1, graph->bndptr);

  nbnd = 0;

  switch (ctrl->objtype) {
    case METIS_OBJTYPE_CUT:
      /* Compute  */
      if (bndtype == BNDTYPE_REFINE) {
        for (i=0; i<nvtxs; i++) {
          if (graph->ckrinfo[i].ed > 0 && graph->ckrinfo[i].ed-graph->ckrinfo[i].id >= 0)
            BNDInsert(nbnd, bndind, bndptr, i);
        }
      }
      else { /* BNDTYPE_BALANCE */
        for (i=0; i<nvtxs; i++) {
          if (graph->ckrinfo[i].ed > 0) 
            BNDInsert(nbnd, bndind, bndptr, i);
        }
      }
      break;

    case METIS_OBJTYPE_VOL:
      /* Compute the boundary */
      if (bndtype == BNDTYPE_REFINE) {
        for (i=0; i<nvtxs; i++) {
          if (graph->vkrinfo[i].gv >= 0)
            BNDInsert(nbnd, bndind, bndptr, i);
        }
      }
      else { /* BNDTYPE_BALANCE */
        for (i=0; i<nvtxs; i++) {
          if (graph->vkrinfo[i].ned > 0) 
            BNDInsert(nbnd, bndind, bndptr, i);
        }
      }
      break;

    default:
      gk_errexit(SIGERR, "Unknown objtype of %d\n", ctrl->objtype);
  }

  graph->nbnd = nbnd;
}



void ComputeKWayVolGains_2(ctrl_t *ctrl, graph_t *graph)
{
  idx_t i, ii, j, k, l, nvtxs, nparts, me, other, pid; 
  idx_t *xadj, *vsize, *adjncy, *adjwgt, *where, 
        *bndind, *bndptr, *ophtable;
  vkrinfo_t *myrinfo, *orinfo;
  vnbr_t *mynbrs, *onbrs;

  WCOREPUSH;

  nparts = ctrl->nparts;

  nvtxs  = graph->nvtxs;
  xadj   = graph->xadj;
  vsize  = graph->vsize;
  adjncy = graph->adjncy;
  adjwgt = graph->adjwgt;

  where  = graph->where;
  bndind = graph->bndind;
  bndptr = iset(nvtxs, -1, graph->bndptr);

  ophtable = iset(nparts, -1, iwspacemalloc(ctrl, nparts));

  graph->minvol = graph->nbnd = 0;
  for (i=0; i<nvtxs; i++) {
    myrinfo     = graph->vkrinfo+i;
    myrinfo->gv = IDX_MIN;

    if (myrinfo->nnbrs > 0) {
      me     = where[i];
      mynbrs = ctrl->vnbrpool + myrinfo->inbr;

      graph->minvol += myrinfo->nnbrs*vsize[i];

      for (j=xadj[i]; j<xadj[i+1]; j++) {
        ii     = adjncy[j];
        other  = where[ii];
        orinfo = graph->vkrinfo+ii;
        onbrs  = ctrl->vnbrpool + orinfo->inbr;

        for (k=0; k<orinfo->nnbrs; k++) 
          ophtable[onbrs[k].pid] = k;
        ophtable[other] = 1;  

        if (me == other) {
         
          for (k=0; k<myrinfo->nnbrs; k++) {
            if (ophtable[mynbrs[k].pid] == -1)
              mynbrs[k].gv -= vsize[ii];
          }
        }
        else {
          ASSERT(ophtable[me] != -1);

          if (onbrs[ophtable[me]].ned == 1) { 
           
            for (k=0; k<myrinfo->nnbrs; k++) {
              if (ophtable[mynbrs[k].pid] != -1) 
                mynbrs[k].gv += vsize[ii];
            }
          }
          else {
            
            for (k=0; k<myrinfo->nnbrs; k++) {
              if (ophtable[mynbrs[k].pid] == -1) 
                mynbrs[k].gv -= vsize[ii];
            }
          }
        }

        
        for (k=0; k<orinfo->nnbrs; k++) 
          ophtable[onbrs[k].pid] = -1;
        ophtable[other] = -1;
      }

      
      for (k=0; k<myrinfo->nnbrs; k++) {
        if (mynbrs[k].gv > myrinfo->gv)
          myrinfo->gv = mynbrs[k].gv;
      }

      
      if (myrinfo->ned > 0 && myrinfo->nid == 0)
        myrinfo->gv += vsize[i];
    }

    if (myrinfo->gv >= 0)
      BNDInsert(graph->nbnd, bndind, bndptr, i);
  }

  WCOREPOP;
}



int IsBalanced_2(ctrl_t *ctrl, graph_t *graph, real_t ffactor)
{
  return 
    (ComputeLoadImbalanceDiff(graph, ctrl->nparts, ctrl->pijbm, ctrl->ubfactors) 
         <= ffactor);
}



void XRK(ctrl_t *ctrl, graph_t *orggraph, graph_t *graph, idx_t *part)
{
  idx_t i, nlevels, contig=ctrl->contig;
  graph_t *ptr;
// for (idx_t i = 0; i < graph->nvtxs; i++) {
//   graph->where[i]=part[i];
// }
// for(idx_t i=0;i<graph->nvtxs;i++){graph->where[i]=part[i];
// }
  IFSET(ctrl->dbglvl, METIS_DBG_TIME, gk_startcputimer(ctrl->UncoarsenTmr));

 
  for (ptr=graph, nlevels=0; ptr==orggraph; ptr=ptr->finer, nlevels++); 

  
  ComputeKWayPartitionParams_2(ctrl, graph);

  
  if (ctrl->minconn) 
    EliminateSubDomainEdges(ctrl, graph);
  
  
  if (contig && FindPartitionInducedComponents(graph, graph->where, NULL, NULL) > ctrl->nparts) { 
    EliminateComponents(ctrl, graph);

    ComputeKWayBoundary_2(ctrl, graph, BNDTYPE_BALANCE);
    GKO(ctrl, graph, 5, 0, OMODE_BALANCE); 

    ComputeKWayBoundary_2(ctrl, graph, BNDTYPE_REFINE);
    GKO(ctrl, graph, ctrl->niter, 0, OMODE_REFINE); 

    ctrl->contig = 0;
  }

 
  for (i=0; ;i++) {
    if (ctrl->minconn && i == nlevels/2) 
      EliminateSubDomainEdges(ctrl, graph);

    IFSET(ctrl->dbglvl, METIS_DBG_TIME, gk_startcputimer(ctrl->RefTmr));

    if (2*i >= nlevels && !IsBalanced_2(ctrl, graph, .02)) {
      ComputeKWayBoundary_2(ctrl, graph, BNDTYPE_BALANCE);
      GKO(ctrl, graph, 1, 0, OMODE_BALANCE); 
      ComputeKWayBoundary_2(ctrl, graph, BNDTYPE_REFINE);
     
    }

    GKO(ctrl, graph, ctrl->niter, 5.0, OMODE_REFINE); 

    IFSET(ctrl->dbglvl, METIS_DBG_TIME, gk_stopcputimer(ctrl->RefTmr));

    
    if (contig && i == nlevels/2) {
      if (FindPartitionInducedComponents(graph, graph->where, NULL, NULL) > ctrl->nparts) {
        EliminateComponents(ctrl, graph);

        if (!IsBalanced_2(ctrl, graph, .02)) {
          ctrl->contig = 1;
          ComputeKWayBoundary_2(ctrl, graph, BNDTYPE_BALANCE);
          GKO(ctrl, graph, 5, 0, OMODE_BALANCE); 
  
          ComputeKWayBoundary_2(ctrl, graph, BNDTYPE_REFINE);
          GKO(ctrl, graph, ctrl->niter, 0, OMODE_REFINE); 
          ctrl->contig = 0;
        }
      }
    }

    if (graph == orggraph)
      break;

    graph = graph->finer;

    graph_ReadFromDisk(ctrl, graph);

    IFSET(ctrl->dbglvl, METIS_DBG_TIME, gk_startcputimer(ctrl->ProjectTmr));
    ASSERT(graph->vwgt != NULL);

    ProjectKWayPartition_2(ctrl, graph);
    IFSET(ctrl->dbglvl, METIS_DBG_TIME, gk_stopcputimer(ctrl->ProjectTmr));
  }

  
  ctrl->contig = contig;
  if (contig && FindPartitionInducedComponents(graph, graph->where, NULL, NULL) > ctrl->nparts) 
    EliminateComponents(ctrl, graph);

  if (!IsBalanced_2(ctrl, graph, 0.0)) {
    ComputeKWayBoundary_2(ctrl, graph, BNDTYPE_BALANCE);
    GKO(ctrl, graph, 10, 0, OMODE_BALANCE); 

    ComputeKWayBoundary_2(ctrl, graph, BNDTYPE_REFINE);
    GKO(ctrl, graph, ctrl->niter, 0, OMODE_REFINE); 
  }
ComputeKWayBoundary_2(ctrl, graph, BNDTYPE_BALANCE);
    GKO(ctrl, graph, 10, 0, OMODE_BALANCE); 

    ComputeKWayBoundary_2(ctrl, graph, BNDTYPE_REFINE);
    GKO(ctrl, graph, ctrl->niter, 0, OMODE_REFINE); 
  if (ctrl->contig) 
    ASSERT(FindPartitionInducedComponents(graph, graph->where, NULL, NULL) == ctrl->nparts);

  IFSET(ctrl->dbglvl, METIS_DBG_TIME, gk_stopcputimer(ctrl->UncoarsenTmr));
}


idx_t adaptiverefine(ctrl_t *ctrl, graph_t *graph, idx_t *part,idx_t *new_part)
{
  idx_t i, j, objval=0, curobj=0, bestobj=0,court=0;
  real_t curbal=0.0, bestbal=0.0;
  graph_t *cgraph;
  int status;
      
  for (idx_t i = 0; i < graph->nvtxs; i++) {
      if (part[i] != 0) {
          court = 1; // 找到非零元素
          break;
      }
  }
  if(court==0)
    {
    for (i=0; i<ctrl->ncuts; i++) {
      
       
        XRK(ctrl, graph, graph, part);
        
        switch (ctrl->objtype) {
          case METIS_OBJTYPE_CUT:
            curobj = graph->mincut;
            break;
          case METIS_OBJTYPE_VOL:
            curobj = graph->minvol;
            break;
          default:
            gk_errexit(SIGERR, "Unknown objtype: %d\n", ctrl->objtype);
        }
        curbal = ComputeLoadImbalanceDiff(graph, ctrl->nparts, ctrl->pijbm, ctrl->ubfactors);
        if (i == 0 
            || (curbal <= 0.0005 && bestobj > curobj)
            || (bestbal > 0.0005 && curbal < bestbal)) {
          // icopy(graph->nvtxs, graph->where, part);
          bestobj = curobj;
          bestbal = curbal;
        }
 for (idx_t i = 0; i < graph->nvtxs; i++) {

       new_part[i]=graph->where[i];   }
        FreeRData(graph);
        if (bestobj == 0)
          break;
        }
    }
  else{
     for (i=0; i<ctrl->ncuts; i++) {
       

        XRK(ctrl, graph, graph, part);
       
        switch (ctrl->objtype) {
          case METIS_OBJTYPE_CUT:
            curobj = graph->mincut;
            break;

          case METIS_OBJTYPE_VOL:
            curobj = graph->minvol;
            break;

          default:
            gk_errexit(SIGERR, "Unknown objtype: %d\n", ctrl->objtype);
        }
        curbal = ComputeLoadImbalanceDiff(graph, ctrl->nparts, ctrl->pijbm, ctrl->ubfactors);
       
    for (idx_t i = 0; i < graph->nvtxs; i++) {
       new_part[i]=graph->where[i];   }
       
        
        FreeRData(graph);
        if (bestobj == 0)
          break;
      }
    }

    
  FreeGraph(&graph);
  return bestobj;
}



void InitKWayPartitioning_2(ctrl_t *ctrl, graph_t *graph)
{
  idx_t i, ntrials, options[METIS_NOPTIONS], curobj=0, bestobj=0;
  idx_t *bestwhere=NULL;
  real_t *ubvec=NULL;
  int status;

  METIS_SetDefaultOptions(options);
  
  options[METIS_OPTION_NITER]     = ctrl->niter;
  options[METIS_OPTION_OBJTYPE]   = METIS_OBJTYPE_CUT;
  options[METIS_OPTION_NO2HOP]    = ctrl->no2hop;
  options[METIS_OPTION_ONDISK]    = ctrl->ondisk;
  options[METIS_OPTION_DROPEDGES] = ctrl->dropedges;
  

  ubvec = rmalloc(graph->ncon, "InitKWayPartitioning_2: ubvec");
  for (i=0; i<graph->ncon; i++) 
    ubvec[i] = (real_t)pow(ctrl->ubfactors[i], 1.0/log(ctrl->nparts));


  switch (ctrl->objtype) {
    case METIS_OBJTYPE_CUT:
    case METIS_OBJTYPE_VOL:
      options[METIS_OPTION_NCUTS] = ctrl->nIparts;
      status = METIS_PartGraphRecursive(&graph->nvtxs, &graph->ncon, 
                   graph->xadj, graph->adjncy, graph->vwgt, graph->vsize, 
                   graph->adjwgt, &ctrl->nparts, ctrl->tpwgts, ubvec, 
                   options, &curobj, graph->where);

      if (status != METIS_OK)
        gk_errexit(SIGERR, "Failed during initial partitioning\n");

      break;

#ifdef XXX 
    case METIS_OBJTYPE_VOL:
      bestwhere = imalloc(graph->nvtxs, "InitKWayPartitioning_2: bestwhere");
      options[METIS_OPTION_NCUTS] = 2;

      ntrials = (ctrl->nIparts+1)/2;
      for (i=0; i<ntrials; i++) {
        status = METIS_PartGraphRecursive(&graph->nvtxs, &graph->ncon, 
                     graph->xadj, graph->adjncy, graph->vwgt, graph->vsize, 
                     graph->adjwgt, &ctrl->nparts, ctrl->tpwgts, ubvec, 
                     options, &curobj, graph->where);
        if (status != METIS_OK)
          gk_errexit(SIGERR, "Failed during initial partitioning\n");

        curobj = ComputeVolume(graph, graph->where);

        if (i == 0 || bestobj > curobj) {
          bestobj = curobj;
          if (i < ntrials-1)
            icopy(graph->nvtxs, graph->where, bestwhere);
        }

        if (bestobj == 0)
          break;
      }
      if (bestobj != curobj)
        icopy(graph->nvtxs, bestwhere, graph->where);

      break;
#endif

    default:
      gk_errexit(SIGERR, "Unknown objtype: %d\n", ctrl->objtype);
  }

  gk_free((void **)&ubvec, &bestwhere, LTERM);

}




 

  
  
 

idx_t GrowMultisection_2(ctrl_t *ctrl, graph_t *graph, idx_t nparts, idx_t *where)
{
  idx_t i, j, k, l, nvtxs, nleft, first, last; 
  idx_t *xadj, *vwgt, *adjncy;
  idx_t *queue;
  idx_t tvwgt, maxpwgt, *pwgts;

  WCOREPUSH;

  nvtxs  = graph->nvtxs;
  xadj   = graph->xadj;
  vwgt   = graph->xadj;
  adjncy = graph->adjncy;

  queue = iwspacemalloc(ctrl, nvtxs);


  
  for (nleft=0, i=0; i<nvtxs; i++) {
    if (xadj[i+1]-xadj[i] > 1) 
      where[nleft++] = i;
  }
  nparts = gk_min(nparts, nleft);
  for (i=0; i<nparts; i++) {
    j = irandInRange(nleft);
    queue[i] = where[j];
    where[j] = --nleft;
  }

  pwgts   = iset(nparts, 0, iwspacemalloc(ctrl, nparts));
  tvwgt   = isum(nvtxs, vwgt, 1);
  maxpwgt = (1.5*tvwgt)/nparts;

  iset(nvtxs, -1, where);
  for (i=0; i<nparts; i++) { 
    where[queue[i]] = i;
    pwgts[i] = vwgt[queue[i]];
  }

  first = 0; 
  last  = nparts;
  nleft = nvtxs-nparts;


  
  while (first < last) { 
    i = queue[first++];
    l = where[i];
    if (pwgts[l] > maxpwgt)
      continue;

    for (j=xadj[i]; j<xadj[i+1]; j++) {
      k = adjncy[j];
      if (where[k] == -1) {
        if (pwgts[l]+vwgt[k] > maxpwgt)
          break;
        pwgts[l] += vwgt[k];
        where[k] = l;
        queue[last++] = k;
        nleft--;
      }
    }
  }
  
  
  if (nleft > 0) { 
    for (i=0; i<nvtxs; i++) {
      if (where[i] == -1)
        where[i] = irandInRange(nparts);
    }
  }

  WCOREPOP;

  return nparts;
}








int newrefine(idx_t *nvtxs, idx_t *ncon, idx_t *xadj, idx_t *adjncy, 
          idx_t *vwgt, idx_t *vsize, idx_t *adjwgt, idx_t *nparts, 
          real_t *tpwgts, real_t *ubvec, idx_t *options, idx_t *objval, 
          idx_t *part, idx_t *new_part)
{
  
  int sigrval=0, renumber=0,court=0;
  graph_t *graph;
  ctrl_t *ctrl;
  idx_t *old_part;
  
  if (!gk_malloc_init()) 
    return METIS_ERROR_MEMORY;

  gk_sigtrap();

  if ((sigrval = gk_sigcatch()) != 0)
    goto SIGTHROW;

  
  ctrl = SetupCtrl(METIS_OP_KMETIS, options, *ncon, *nparts, tpwgts, ubvec);
  if (!ctrl) {
    gk_siguntrap();
    return METIS_ERROR_INPUT;
  }

  
  if (ctrl->numflag == 1) {
    Change2CNumbering(*nvtxs, xadj, adjncy);
    renumber = 1;
  }

  
  graph = SetupGraph(ctrl, *nvtxs, *ncon, xadj, adjncy, vwgt, vsize, adjwgt);

  
  SetupKWayBalMultipliers(ctrl, graph);

  
  ctrl->CoarsenTo = gk_max((*nvtxs)/(40*gk_log2(*nparts)), 30*(*nparts));
  ctrl->nIparts   = (ctrl->nIparts != -1 ? ctrl->nIparts : (ctrl->CoarsenTo == 30*(*nparts) ? 4 : 5));

  
  if (ctrl->contig && !IsConnected(graph, 0)) 
    gk_errexit(SIGERR, "METIS Error: A contiguous partition is requested for a non-contiguous input graph.\n");
    
    
  AllocateWorkSpace(ctrl, graph);

  
  IFSET(ctrl->dbglvl, METIS_DBG_TIME, InitTimers(ctrl));
  IFSET(ctrl->dbglvl, METIS_DBG_TIME, gk_startcputimer(ctrl->TotalTmr));

  for (idx_t i = 0; i < graph->nvtxs; i++) {
      if (part[i] != 0) {
          court = 1; 
          break;
      }
  }

  if(court==0)
    {iset(*nvtxs, 0, part);
    for (idx_t i=0; i<ctrl->ncuts; i++) {
      
        IFSET(ctrl->dbglvl, METIS_DBG_TIME, gk_startcputimer(ctrl->InitPartTmr));
        AllocateKWayPartitionMemory_2(ctrl, graph);
        FreeWorkSpace(ctrl);
        InitKWayPartitioning_2(ctrl, graph);
        AllocateWorkSpace(ctrl, graph);
        AllocateRefinementWorkSpace(ctrl, graph->nedges, graph->nedges);
        IFSET(ctrl->dbglvl, METIS_DBG_TIME, gk_stopcputimer(ctrl->InitPartTmr));
        IFSET(ctrl->dbglvl, METIS_DBG_IPART, 
            printf("Initial %ld"PRIDX"-way partitioning cut: %ln"PRIDX"\n", ctrl->nparts, objval));}}
  else{
    graph->where = (idx_t *)malloc(graph->nvtxs * sizeof(idx_t));
      
        
      for (idx_t i=0; i<ctrl->ncuts; i++) {
      
        IFSET(ctrl->dbglvl, METIS_DBG_TIME, gk_startcputimer(ctrl->InitPartTmr));
        AllocateKWayPartitionMemory_2(ctrl, graph);
        FreeWorkSpace(ctrl);
       
        AllocateWorkSpace(ctrl, graph);
        AllocateRefinementWorkSpace(ctrl, graph->nedges, graph->nedges);
        IFSET(ctrl->dbglvl, METIS_DBG_TIME, gk_stopcputimer(ctrl->InitPartTmr));
        IFSET(ctrl->dbglvl, METIS_DBG_IPART, 
            printf("Initial %ld"PRIDX"-way partitioning cut: %ln"PRIDX"\n", ctrl->nparts, objval));}
            for (idx_t i = 0; i < graph->nvtxs; i++) {
          graph->where[i] = part[i];
      }
           
    }
    
    
  *objval = (*nparts == 1 ? 0 : adaptiverefine(ctrl, graph, part,new_part));
  IFSET(ctrl->dbglvl, METIS_DBG_TIME, gk_stopcputimer(ctrl->TotalTmr));
  IFSET(ctrl->dbglvl, METIS_DBG_TIME, PrintTimers(ctrl));
    FreeCtrl(&ctrl);

SIGTHROW:
  if (renumber)
    Change2FNumbering(*nvtxs, xadj, adjncy, part);

  gk_siguntrap();
  gk_malloc_cleanup(0);

  return metis_rcode(sigrval);
}
