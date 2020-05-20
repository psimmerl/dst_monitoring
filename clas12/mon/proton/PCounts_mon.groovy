package clas12.mon.proton

import org.jlab.detector.base.DetectorType
import org.jlab.clas.physics.LorentzVector
import org.jlab.groot.data.H1F
import org.jlab.groot.data.H2F
import java.util.concurrent.ConcurrentHashMap


class PCounts_mon {
  def hists = new ConcurrentHashMap()
  def entry = new ConcurrentHashMap()

  def processEvent(event) {
    if(event.hasBank("RUN::config")) {
      def cnfb = event.getBank("RUN::config")
      def ts = cnfb.getLong("timestamp", 0)
      def evn = cnfb.getInt("event", 0)

      if(event.hasBank("RUN::scaler")) {
        def scaler = event.getBank("RUN::scaler")
        def fcg = scaler.getFloat('fcupgated',0)
        entry[evn] = [ts, 0, fcg]
      }

     if(event.hasBank("REC::Particle") && event.hasBank("REC::Calorimeter")) {
        def partb = event.getBank("REC::Particle")
        def calb = event.getBank("REC::Calorimeter")
        def ipro = (0..<partb.rows()).find{partb.getInt('pid',it)==2212 && partb.getShort("status",it)>=400}# also found partb.getShort("status",it)<4000 and no partb.getShort dependence
        def psec = (0..<calb.rows()).findResult{
          (calb.getShort('pindex',it).toInteger() == ipro &&
          calb.getByte('detector',it).toInteger() == DetectorType.ECAL.getDetectorId()) ? calb.getByte('sector',it) : null
        }

        if(ipro!=null && psec!=null)
          entry[evn+0.1] = [ts, 1, psec]
      }
    }
  }


  def finish() {
    def tline = entry.sort{it.key}.collect{it.value}

    def data = [[fcup: [], npro: [0]*6]]
    tline.each{ts,id,val->
      if(id==0) {
        def fcg = val
        data[-1].fcup.add([ts,fcg])

        data[-1].dt = data.last().fcup.with{last()[0]-first()[0]}
        data[-1].dq = data.last().fcup.with{last()[1]-first()[1]}
        data[-1].norm = data[-1].with{npro*.div(dq)}

        if(data[-1].fcup.with{last()[0]-first()[0]>4e8}) {
          data.add([fcup: [[ts,fcg]], npro: [0]*6])
        }
      } else if(data[-1].fcup) {
        def isec = val-1
        data[-1].npro[isec]++
      }
    }

    int maxnpro = data.collect{it.npro.max()}.max().toInteger()+10
    def maxnorm = data.dropRight(1).collect{it.norm.max()}.max()*1.05
    def hsnpro = (1..6).collect{new H2F("h2npro_s${it}", "number of protons in sec $it between FC readings;number of protons", maxnpro,0,maxnpro,200,0,60)}
    def hsnorm0 = (1..6).collect{new H2F("full/h2pronorm0_s${it}", "normalized number of protons in sec $it;normalized number of protons", 200,0,maxnorm,200,0,60)}
    def hsnorm1 = (1..6).collect{new H2F("fixed/h2pronorm1_s${it}", "normalized number of protons in sec $it;normalized number of protons", 200,0,5,200,0,60)}

    maxnpro = data.collect{it.npro.sum()}.max().toInteger()+10
    maxnorm = data.dropRight(1).collect{it.norm.sum()}.max()*1.05
    def h0npro = new H2F("h2npro", "number of protons in all sectors between FC readings;number of protons", maxnpro,0,maxnpro,200,0,60)
    def h0norm0 = new H2F("full/h2pronorm0", "normalized number of protons in all sectors;normalized number of protons", 200,0,maxnorm,200,0,60)
    def h0norm1 = new H2F("fixed/h2pronorm1", "normalized number of protons in all sectors;normalized number of protons", 200,0,25,200,0,60)

    data.dropRight(1).each{
      def curr = it.dq/it.dt/4*1e9
      (0..<6).each{i->hsnpro[i].fill(it.npro[i], curr)}
      (0..<6).each{i->hsnorm0[i].fill(it.norm[i], curr)}
      (0..<6).each{i->hsnorm1[i].fill(it.norm[i], curr)}
      h0npro.fill(it.npro.sum(), curr)
      h0norm0.fill(it.norm.sum(), curr)
      h0norm1.fill(it.norm.sum(), curr)
    }

    [hsnpro,hsnorm0,hsnorm1].each{it.each{hists[it.getName().replace("h2","h")] = it.projectionX(it.getName().replace("h2","h"), it.getTitle())}}
    [h0npro,h0norm0,h0norm1].each{hists[it.getName().replace('h2','h')] = it.projectionX(it.getName().replace("h2","h"), it.getTitle())}

    [hsnpro,hsnorm0,hsnorm1].each{it.each{hists[it.getName()] = it}}
    [h0npro,h0norm0,h0norm1].each{hists[it.getName()] = it}
  }
}
