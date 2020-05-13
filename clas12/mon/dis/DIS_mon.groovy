package clas12.mon.dis

import org.jlab.clas.physics.LorentzVector
import org.jlab.groot.data.H1F
import org.jlab.groot.data.H2F
import java.util.concurrent.ConcurrentHashMap


class DIS_mon {
  def hists = new ConcurrentHashMap()

  def beam = LorentzVector.withPID(11,0,0,10.6)
  def target = LorentzVector.withPID(2212,0,0,0)

  def hw0 = {new H1F("$it","$it",200,0.5,2.5)}
  def hww = {new H1F("$it","$it",200,0,5)}
  def hthepe = {new H2F("$it","$it",200,0,11,200,0,40)}

  def E0 = 10.6

  def banknames = ['REC::Particle','REC::Calorimeter']
  def processEvent(event) {
    if(banknames.every{event.hasBank(it)}) {
      def (partb,ecb) = banknames.collect{event.getBank(it)}

      (0..<partb.rows()).find{partb.getInt('pid',it)==11 && partb.getShort("status",it)<0}
        ?.with{iele->
          def ele = LorentzVector.withPID(11,*['px','py','pz'].collect{partb.getFloat(it,iele)})
          def chi2pid = partb.getFloat('chi2pid',iele)

          def esec = (0..<ecb.rows()).find{ecb.getShort('pindex',it)==0}.with{ecb.getByte('sector',it)}

          def ex = beam+target-ele

          def yy = (E0-ele.e())/E0
          def wname = ex.mass()<2 ? 'w.lt.2' : 'w.gt.2'

          if(ex.mass()>1.1)
          [
           ['y.lt.1', yy<1],
           ['y.lt.0.85', yy<0.85],
           ['y.lt.0.8', yy<0.8],
           ['y.lt.0.7', yy<0.7],
           ['y.lt.0.75', yy<0.75],
           ['chi2pid.lt.3', chi2pid<3],
           ['y.lt.1/chi2pid.lt.3', chi2pid<3 && yy<1],
           ['y.lt.0.85/chi2pid.lt.3', chi2pid<3 && yy<0.85],
           ['y.lt.0.8/chi2pid.lt.3', chi2pid<3 && yy<0.8],
           ['y.lt.0.7/chi2pid.lt.3', chi2pid<3 && yy<0.7],
           ['y.lt.0.75/chi2pid.lt.3', chi2pid<3 && yy<0.75],
          ].findAll{it[1]}.each{name,cut->
            hists.computeIfAbsent("$name/hW0",hw0).fill(ex.mass())
            hists.computeIfAbsent("$name/hW_s$esec",hw0).fill(ex.mass())
            hists.computeIfAbsent("$name/hW",hww).fill(ex.mass())
            hists.computeIfAbsent("$wname/$name/hW",hww).fill(ex.mass())
            hists.computeIfAbsent("$wname/$name/hW_$esec",hww).fill(ex.mass())
            hists.computeIfAbsent("$wname/$name/hthepe",hthepe).fill(ele.e(), Math.toDegrees(ele.theta()))
          }
        }
    }
  }
}
