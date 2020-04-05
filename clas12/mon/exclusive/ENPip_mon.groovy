package clas12.mon.exclusive

import org.jlab.clas.physics.LorentzVector
import org.jlab.groot.data.H1F
import org.jlab.groot.data.H2F
import java.util.concurrent.ConcurrentHashMap


class ENPip_mon {
  def hists = new ConcurrentHashMap()

  def beam = LorentzVector.withPID(11,0,0,10.6)
  def target = LorentzVector.withPID(2212,0,0,0)

  def hmmneu = {new H1F("$it","$it",200,0,2)}
  def hmm2pip = {new H1F("$it","$it",200,-1,1)}

  def banknames = ['REC::Particle','REC::Calorimeter']
  def processEvent(event) {
    if(banknames.every{event.hasBank(it)}) {
      def (partb,ecb) = banknames.collect{event.getBank(it)}

      (0..<partb.rows()).findAll{partb.getInt('pid',it)==11 && partb.getShort("status",it)<0}
        .collectMany{iele->
          (0..<partb.rows()).findAll{partb.getInt('pid',it)==2112}.collect{ineu->[iele,ineu]}
        }.collectMany{iele,ineu->
          (0..<partb.rows()).findAll{partb.getInt('pid',it)==211 && partb.getShort('status',it)<4000}.collect{ipip->[iele,ineu,ipip]}
        }.each{iele,ineu,ipip->
          def ele = LorentzVector.withPID(11,*['px','py','pz'].collect{partb.getFloat(it,iele)})
          def neu = LorentzVector.withPID(2112,*['px','py','pz'].collect{partb.getFloat(it,ineu)})
          def pip = LorentzVector.withPID(211,*['px','py','pz'].collect{partb.getFloat(it,ipip)})

          def esec = (0..<ecb.rows()).find{ecb.getShort('pindex',it)==0}.with{ecb.getByte('sector',it)}

          def misspip = beam+target-ele-neu
          def missneu = beam+target-ele-pip

          def q2 = -(beam-ele).mass2()

          hists.computeIfAbsent("mm2_pip_$esec",hmm2pip).fill(misspip.mass2())
          hists.computeIfAbsent("mm_neu_$esec",hmmneu).fill(missneu.mass())
        }
    }
  }
}
