package clas12.mon.exclusive

import org.jlab.clas.physics.LorentzVector
import org.jlab.groot.data.H1F
import org.jlab.groot.data.H2F
import java.util.concurrent.ConcurrentHashMap


class EPipNspectator_mon {
  def hists = new ConcurrentHashMap()

  def beam = LorentzVector.withPID(11,0,0,10.6)
  def target = LorentzVector.withPID(2212,0,0,0)

  def hmmneu = {new H1F("$it","$it",200,0.5,2)}
  def hepth = {new H2F("$it","$it",200,0,10,200,0,40)}

  def banknames = ['REC::Particle','REC::Calorimeter']
  def processEvent(event) {
    if(banknames.every{event.hasBank(it)}) {
      def (partb,ecb) = banknames.collect{event.getBank(it)}

      (0..<partb.rows()).findAll{partb.getInt('pid',it)==11 && partb.getShort("status",it)<0}
        .collectMany{iele->
          (0..<partb.rows()).findAll{partb.getInt('pid',it)==211 && partb.getFloat("chi2pid",it)<2}.collect{ipip->[iele,ipip]}
        }.each{iele,ipip->
          def ele = LorentzVector.withPID(11,*['px','py','pz'].collect{partb.getFloat(it,iele)})
          def pip = LorentzVector.withPID(211,*['px','py','pz'].collect{partb.getFloat(it,ipip)})

          def yy = (beam.e()-ele.e())/beam.e()

          def pdet = partb.getShort('status',ipip)<4000 ? 'FD' : 'CD'

          def ex = beam+target-ele
          def epipx = beam+target-ele-pip

          def esec = (0..<ecb.rows()).find{ecb.getShort('pindex',it)==iele}.with{ecb.getByte('sector',it)}

          def wname = ex.mass()<2 ? 'w.lt.2' : 'w.gt.2'
//          if(epipx.mass()<1.04 && epipx.mass()>0.8) {
            hists.computeIfAbsent("$wname/hepth:$pdet",hepth).fill(ele.p(), Math.toDegrees(ele.theta()))
            hists.computeIfAbsent("mm_neu:$pdet:s$esec",hmmneu).fill(epipx.mass())
//          }
        }
    }
  }
}
