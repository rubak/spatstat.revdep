#include "arcs.h"

std::vector<std::vector<double> > morphoArcs(Graph *graph)
{
  if(graph->dbg)Rprintf("mArcs[");
  std::vector<double> value;
  std::vector<std::vector<double> > arcs;
  std::vector<std::vector<double> > arcsi, arcsi2;
  std::vector<double> *parc;
  int i,j,k;
  double angleij;
  double ang0;
  double R = graph->par; // the circles are assumed to be of radius R/2
  double Rcheck = R/2;
  if(R>0) {
    // if R>0
    for(i=0; i<graph->pp->size(); i++)// compute arcs that overlap circle b(x_i, R/2) 
    {//if(graph->pp->included(&i, &Rcheck)){// border correction
      arcsi.clear();
      arcsi2.clear();
      parc = new std::vector<double>;
      parc->resize(3);
      parc->at(0)=(double) i;
      for(j=0; j<(int)graph->nodelist.at(i).size(); j++) {
        k = graph->nodelist.at(i).at(j)-1;
        {  //if(graph->pp->included(&k, &Rcheck)){// border correction not needed here
          angleij = atan2(graph->pp->getY(&k)-graph->pp->getY(&i), graph->pp->getX(&k)-graph->pp->getX(&i));
          if(angleij<0) angleij=angleij+2*PI; // map to [0,2PI]
          parc->at(2)=acos( graph->pp->getDistance(&i, &k) /R); // not 2R but 2R/2. see above.
          parc->at(1)= -1.0*(parc->at(2))+angleij; // clockwise order of start-end
          if(parc->at(1)<0) parc->at(1)= 2*PI+parc->at(1);
          parc->at(2)=parc->at(2)+angleij;
          if(parc->at(2)>2*PI)parc->at(2)=parc->at(2)-2*PI;
          if(parc->at(2)<parc->at(1)){ // if we go over 2PI, split to two
            ang0 = parc->at(2);
            parc->at(2)=2*PI;
            arcsi.push_back(*parc);
            parc->at(1)=0;
            parc->at(2)=ang0;
          }
          arcsi.push_back(*parc);
        }// border correction
      }
      // first arc is the full circle
      parc->at(1) = 0;
      parc->at(2) = 2*PI;
      arcsi2.clear();
      arcsi2.push_back(*parc);
      delete parc;
      if((int)arcsi.size()>0) { // split the cirle to sub arcs
        for(j=0; j<(int)arcsi.size(); j++) {
          morphoArcsMinus(&arcsi2, arcsi.at(j));
        }
      }
      for(j=0; j<(int)arcsi2.size();j++) {
        arcs.push_back(arcsi2.at(j));
      }
    } // border correction
    // end of R>0
  }
  if(graph->dbg)Rprintf("]");
  return arcs;
}

/*********************************************************************************
 * * setminus an arc from a given set of arcs. replaces first argument
 * * *********************************************************************************/
void  morphoArcsMinus(std::vector<std::vector<double> > *arcs, std::vector<double> arc) {
  int i;
  std::vector<std::vector<double> > arcs2;
  std::vector<double> *parc;
  parc = new std::vector<double >; parc->resize(3);
  parc->at(0)=arc.at(0);
  for(i=0; i < (int)arcs->size(); i++) { // split each arc
    if( (arc.at(1) <= arcs->at(i).at(2)) && (arc.at(2) >=  arcs->at(i).at(1)) ) { // if overlap occurs
      if(arc.at(2) >= arcs->at(i).at(2)) { // cut from the end
        if(arc.at(1) > arcs->at(i).at(1)) { // not all is cut
          parc->at(1)=arcs->at(i).at(1);
          parc->at(2)=arc.at(1);
          arcs2.push_back(*parc);
        }
      }
      else if(arc.at(1)<=arcs->at(i).at(1)){ // not cut from the end, but from beginning
        parc->at(1)=arc.at(2);
        parc->at(2)=arcs->at(i).at(2);
        arcs2.push_back(*parc);
      }
      else { // cut from the middle
        parc->at(1)=arcs->at(i).at(1);// first arc
        parc->at(2)=arc.at(1);
        arcs2.push_back(*parc);
        parc->at(1)=arc.at(2);// second arc
        parc->at(2)=arcs->at(i).at(2);
        arcs2.push_back(*parc);
      }
    }
    else { // no overlap
      arcs2.push_back(arcs->at(i));
    }
  }
  delete parc;
  arcs->swap(arcs2);
}


