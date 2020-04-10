package clas12.mon.ft

import org.jlab.detector.base.DetectorType
import org.jlab.clas.physics.LorentzVector
import org.jlab.groot.data.GraphErrors
import org.jlab.groot.data.H1F
import org.jlab.groot.data.H2F
import java.util.concurrent.ConcurrentHashMap


class FT_mon {
  def hists = new ConcurrentHashMap()
  def entry = new ConcurrentHashMap()

  def processEvent(event) {
    if(event.hasBank("RUN::config")) {
      def cnfb = event.getBank("RUN::config")
      def ts = cnfb.getLong("timestamp", 0)
      def evn = cnfb.getInt("event", 0)

      def trg = cnfb.getLong("trigger", 0)
      int trg24 = (trg>>24)&0x1, trg25=(trg>>25)&0x1
      if(trg24 || trg25)
        entry[evn+0.1] = [ts, 1, [trg24, trg25]]

      if(event.hasBank("RUN::scaler")) {
        def scaler = event.getBank("RUN::scaler")
        def fcg = scaler.getFloat('fcupgated',0)
        entry[evn] = [ts, 0, fcg]
      }

      if(event.hasBank("REC::Particle")) {
        def partb = event.getBank("REC::Particle")
        def ifts = (0..<partb.rows()).findAll{partb.getShort("status",it).with{it>=1000 && it<2000} && Math.sqrt(['px','py','pz'].sum{comp->partb.getFloat(comp,it)**2})>0.3}
        int nneuts = ifts.findAll{partb.getByte('charge',it)==0}.size()
        int ncharged = ifts.size() - nneuts

        entry[evn+0.2] = [ts, 2, [nneuts, ncharged]]
      }
    }
  }


  def finish() {
    def data = []
    def tline = entry.sort{it.key}.each{
      def (ts,id,val) = it.value

      if(id==0) {
        def fcg = val

        if(data) {
          data[-1].dt = ts - data[-1].ts
          data[-1].tt = ts - data[0].ts
          data[-1].dq = fcg - data[-1].fc
          data[-1].norm = data[-1].with{ntrg*.div(dq)}
        }

        data.add([fc:fcg, ts:ts, ntrg:[0,0], nqqs: [0,0]])
      } else if(id==1 && data) {
        data[-1].ntrg[0] += val[0]
        data[-1].ntrg[1] += val[1]
      } else if(id==2 && data) {
        data[-1].nqqs[0] += val[0]
        data[-1].nqqs[1] += val[1]
      }
    }

    def grnqqs = [new GraphErrors("grnqqs0"), new GraphErrors("grnqqs1")]
    def grnorm = [new GraphErrors("grnorm24"), new GraphErrors("grnorm25")]
    data.dropRight(1).groupBy{it.tt/4/1e9/30}.each{tt,evlist->
      grnorm[0].addPoint(tt*30, evlist.isEmpty()?0:evlist.sum{it.norm[0]}/evlist.size(), 0,0)
      grnorm[1].addPoint(tt*30, evlist.isEmpty()?0:evlist.sum{it.norm[1]}/evlist.size(), 0,0)
      grnqqs[0].addPoint(tt*30, evlist.isEmpty()?0:evlist.sum{it.nqqs[0]}/evlist.size(), 0,0)
      grnqqs[1].addPoint(tt*30, evlist.isEmpty()?0:evlist.sum{it.nqqs[1]}/evlist.size(), 0,0)
    }

    def maxnqqs = (0..1).collect{i->data.max{it.nqqs[i]}.nqqs[i].toInteger()+5}.collect{[it,0,it]}
    def maxntrg = (0..1).collect{i->data.max{it.ntrg[i]}.ntrg[i].toInteger()+5}.collect{[it,0,it]}
    def maxnorm = (0..1).collect{i->(data.dropRight(1).max{it.norm[i]}.norm[i].toInteger()+5).with{[it*(i?:10),0,it]}}
    def fxdmax = [24:15, 25:200]

    def hsnqqs = [0,1].collect{new H2F("h2nqqs_${it}", "number of ${it?'charged':'neutrals'} between FC readings;number of particles", *maxnqqs[it],200,0,60)}
    def hsntrgs = [24,25].collect{new H2F("h2ntrg_b${it}", "number of ${it}th trigger bits between FC readings;number of trigger bits", *maxntrg[it-24],200,0,60)}
    def hsnorm0 = [24,25].collect{new H2F("full/h2norm0_b${it}", "normalized number of ${it}th trigger bits [full range];normalized number of trigger bits",*maxnorm[it-24],200,0,60)}
    def hsnorm1 = [24,25].collect{new H2F("fixed/h2norm1_b${it}", "normalized number of ${it}th trigger bits [fixed range];normalized number of trigger bits", 200,0,fxdmax[it],200,0,60)}

    data.dropRight(1).each{
      def curr = it.dq/it.dt/4*1e9
      it.nqqs.eachWithIndex{nn,i->hsnqqs[i].fill(nn,curr)}
      it.ntrg.eachWithIndex{nn,i->hsntrgs[i].fill(nn,curr)}
      it.norm.eachWithIndex{nn,i->
        hsnorm0[i].fill(nn,curr)
        hsnorm1[i].fill(nn,curr)
      }
    }

    [hsnqqs, hsntrgs,hsnorm0,hsnorm1].each{it.each{hists[it.getName().replace("h2","h")] = it.projectionX(it.getName().replace("h2","h"), it.getTitle())}}
    [hsnqqs, hsntrgs,hsnorm0,hsnorm1,grnorm,grnqqs].each{it.each{hists[it.getName()] = it}}
  }
}
