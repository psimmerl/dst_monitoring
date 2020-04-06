package clas12.mon.ft

import org.jlab.detector.base.DetectorType
import org.jlab.clas.physics.LorentzVector
import org.jlab.groot.data.H1F
import org.jlab.groot.data.H2F
import java.util.concurrent.ConcurrentHashMap


class FT_mon {
  def hists = new ConcurrentHashMap()
  def fcmap = new ConcurrentHashMap()
  def trgmap = new ConcurrentHashMap()

  def processEvent(event) {
    if(event.hasBank("RUN::config")) {
      def cnfb = event.getBank("RUN::config")
      def ts = cnfb.getLong("timestamp", 0)
      def trg = cnfb.getLong("trigger", 0)
      int trg24 = (trg>>24)&0x1, trg25=(trg>>25)&0x1
      if(trg24 || trg25) trgmap[ts] = [trg24, trg25]

      if(event.hasBank("RUN::scaler")) {
        def scaler = event.getBank("RUN::scaler")
        def fcg = scaler.getFloat('fcupgated',0)
        fcmap[ts] = fcg
      }

/*
     if(event.hasBank("REC::Particle") && event.hasBank("REC::Calorimeter")) {
        def partb = event.getBank("REC::Particle")
        def calb = event.getBank("REC::Calorimeter")
        def iele = (0..<partb.rows()).find{partb.getInt('pid',it)==11 && partb.getShort("status",it)<0}
        def esec = (0..<calb.rows()).findResult{
          (calb.getShort('pindex',it).toInteger() == iele &&
          calb.getByte('detector',it).toInteger() == DetectorType.ECAL.getDetectorId()) ? calb.getByte('sector',it) : null
        }

        if(iele!=null && esec!=null) elentry[ts] = esec
      }
*/
    }
  }


  def finish() {
    def tline = fcmap.collect{ts,fcg->[ts,0,fcg]}
    tline.addAll(trgmap.collect{ts,trgs->[ts,1,trgs]})
    tline.sort{a,b->a[0]<=>b[0]?:a[1]<=>b[1]}

    def data = []
    tline.each{ts,id,val->
      if(id==0) {
        def fcg = val
        if(data) {
          def (ts0, fc0, trgs) = data.last()
          def dt = ts - ts0
          def dq = fcg - fc0
          data[-1].addAll([dt, dq, trgs*.div(dq)])
        }
        data.add([ts,fcg,[0,0]])
      } else if(data) {
        data[-1][2][0] += val[0]
        data[-1][2][1] += val[1]
      }
    }

    int maxntrg = data.collect{it[2].max()}.max().toInteger()+5
    def maxnorm = data.dropRight(1).collect{it[5].max()}.max()*1.1
    def fxdmax = [24:15, 25:200]
    def hsntrgs = [24,25].collect{new H2F("h2ntrg_b${it}", "number of ${it}th trigger bits between FC readings;number of trigger bits", maxntrg,0,maxntrg,200,0,60)}
    def hsnorm0 = [24,25].collect{new H2F("full/h2norm_b${it}", "normalized number of ${it}th trigger bits [full range];normalized number of trigger bits", 200,0,maxnorm,200,0,60)}
    def hsnorm1 = [24,25].collect{new H2F("fixed/h2norm_b${it}", "normalized number of ${it}th trigger bits [fixed range];normalized number of trigger bits", 200,0,fxdmax[it],200,0,60)}

    data.dropRight(1).each{ts,fcg,ntrg,dts,dfc,norm->
      def curr = dfc/dts/4*1e9
      ntrg.eachWithIndex{nn,i->hsntrgs[i].fill(nn,curr)}
      norm.eachWithIndex{nn,i->
        hsnorm0[i].fill(nn,curr)
        hsnorm1[i].fill(nn,curr)
      }
    }

    [hsntrgs,hsnorm0,hsnorm1].each{it.each{hists[it.getName().replace("h2","h")] = it.projectionX(it.getName().replace("h2","h"), h2.getTitle())}}
    [hsntrgs,hsnorm0,hsnorm1].each{it.each{hists[it.getName()] = it}}
  }
}
