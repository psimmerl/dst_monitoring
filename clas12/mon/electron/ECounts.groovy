package clas12.mon.electron

import org.jlab.detector.base.DetectorType
import org.jlab.clas.physics.LorentzVector
import org.jlab.groot.data.H1F
import org.jlab.groot.data.H2F
import java.util.concurrent.ConcurrentHashMap


class ECounts {
  def hists = new ConcurrentHashMap()
  def fcentry = new ConcurrentHashMap()
  def elentry = new ConcurrentHashMap()

  def processEvent(event) {
    if(event.hasBank("RUN::config")) {
      def cnfb = event.getBank("RUN::config")
      def ts = cnfb.getLong("timestamp", 0)

      if(event.hasBank("RUN::scaler")) {
        def scaler = event.getBank("RUN::scaler")
        def fcg = scaler.getFloat('fcupgated',0)
        fcentry[ts] = fcg
      }

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
    }
  }


  def finish() {
    def tline = fcentry.collect{ts,fcg->[ts,0,fcg]}
    tline.addAll(elentry.collect{ts,esec->[ts,1,esec]})
    tline.sort{a,b->a[0]<=>b[0]?:a[1]<=>b[1]}

    def data = []
    tline.each{ts,id,val->
      if(id==0) {
        def fcg = val
        if(data) {
          def (ts0, fc0, nele) = data.last()
          def dt = ts - ts0
          def dq = fcg - fc0
          data[-1].addAll([dt, dq, nele.collect{it/dq}])
        }
        data.add([ts,fcg,[0]*6])
      } else if(data) {
        def isec = val-1
        data[-1][2][isec]++
      }
    }

    int maxnele = Math.max(20, data.collect{it[2].max()}.max().toInteger()+5)
    def maxnorm = data.dropRight(1).collect{it[5].max()}.max()
    def hsnele = (1..6).collect{new H2F("h2nele_s${it}", "number of electrons in sec $it between FC readings;number of electrons", maxnele,0,maxnele,200,0,60)}
    def hsnorm0 = (1..6).collect{new H2F("full/h2enorm_s${it}", "normalized number of electrons in sec $it;normalized number of electrons", 200,0,maxnorm,200,0,60)}
    def hsnorm1 = (1..6).collect{new H2F("fixed/h2enorm_s${it}", "normalized number of electrons in sec $it;normalized number of electrons", 200,0,15,200,0,60)}

    maxnele = Math.max(50, data.collect{it[2].sum()}.max().toInteger()+5)
    maxnorm = data.dropRight(1).collect{it[5].sum()}.max()
    def h0nele = new H2F("h2nele", "number of electrons in all sectors between FC readings;number of electrons", maxnele,0,maxnele,200,0,60)
    def h0norm0 = new H2F("full/h2enorm", "normalized number of electrons in all sectors;normalized number of electrons", 200,0,maxnorm,200,0,60)
    def h0norm1 = new H2F("fixed/h2enorm", "normalized number of electrons in all sectors;normalized number of electrons", 200,0,40,200,0,60)

    data.dropRight(1).each{ts,fcg,nele,dts,dfc,norm->
      def curr = dfc/dts/4*1e9
      (0..<6).each{hsnele[it].fill(nele[it],curr)}
      (0..<6).each{hsnorm0[it].fill(norm[it],curr)}
      (0..<6).each{hsnorm1[it].fill(norm[it],curr)}
      h0nele.fill(nele.sum(),curr)
      h0norm0.fill(norm.sum(),curr)
      h0norm1.fill(norm.sum(),curr)
    }

    [hsnele,hsnorm0,hsnorm1].each{it.each{hists[it.getName().replace("h2","h")] = it.projectionX(it.getName().replace("h2","h"))}}
    [h0nele,h0norm0,h0norm1].each{hists[it.getName().replace('h2','h')] = it.projectionX(it.getName().replace("h2","h"))}

    [hsnele,hsnorm0,hsnorm1].each{it.each{hists[it.getName()] = it}}
    [h0nele,h0norm0,h0norm1].each{hists[it.getName()] = it}
  }
}
