package clas12.mon.exclusive

import org.jlab.clas.physics.LorentzVector
import org.jlab.groot.data.H1F
import org.jlab.groot.data.H2F
import java.util.concurrent.ConcurrentHashMap


class EPPi0_mon {
  def hists = new ConcurrentHashMap()

  def beam = LorentzVector.withPID(11,0,0,10.6)
  def target = LorentzVector.withPID(2212,0,0,0)

  def hmm2epx = {new H1F("$it","$it",200,-1,2)}
  def hggm = {new H1F("$it","$it",200,0.05,0.2)}

  def banknames = ['REC::Particle','REC::Calorimeter']
  def processEvent(event) {
    if(banknames.every{event.hasBank(it)}) {
      def (partb,ecb) = banknames.collect{event.getBank(it)}

      def igs = (0..<partb.rows()).findAll{partb.getInt('pid',it)==22 && partb.getShort('status',it)>=2000}
      def ggs = [igs,igs].combinations().findAll{it[1]>it[0]}
        .findResults{ig1,ig2->
          def g1 = LorentzVector.withPID(22,*['px','py','pz'].collect{partb.getFloat(it,ig1)})
          def g2 = LorentzVector.withPID(22,*['px','py','pz'].collect{partb.getFloat(it,ig2)})
          def ggmass = (g1+g2).mass()
          hists.computeIfAbsent("hggm",hggm).fill(ggmass)
          return (ggmass<0.2 && ggmass>0.05) ? [ig1,ig2,ggmass] : null
        }

      def eps = (0..<partb.rows()).findAll{partb.getInt('pid',it)==11 && partb.getShort("status",it)<0}
        .collectMany{iele->
          (0..<partb.rows()).findAll{partb.getInt('pid',it)==2212 && partb.getShort('status',it)>=4000}.collect{ipro->[iele,ipro]}
        }.collect{iele,ipro->
          def ele = LorentzVector.withPID(11,*['px','py','pz'].collect{partb.getFloat(it,iele)})
          def pro = LorentzVector.withPID(2212,*['px','py','pz'].collect{partb.getFloat(it,ipro)})
          def epx = beam+target-ele-pro

          hists.computeIfAbsent("hepx",hmm2epx).fill(epx.mass2())

          def pi0s = ggs.findAll{ig1,ig2,ggmass->
            def g1 = LorentzVector.withPID(22,*['px','py','pz'].collect{partb.getFloat(it,ig1)})
            def g2 = LorentzVector.withPID(22,*['px','py','pz'].collect{partb.getFloat(it,ig2)})
            def gg = g1+g2
            def epggx = epx-gg
            return ele.vect().theta(g1.vect())>8 && ele.vect().theta(g2.vect())>8 &&
              epx.vect().theta(gg.vect())<2 && epx.mass2()<1 && epggx.e()<1
          }
          return [iele,ipro,pi0s]
        }

      eps.findAll{it[2]}.each{iele,ipro,pi0s->
        def ele = LorentzVector.withPID(11,*['px','py','pz'].collect{partb.getFloat(it,iele)})
        def pro = LorentzVector.withPID(2212,*['px','py','pz'].collect{partb.getFloat(it,ipro)})
        def epx = beam+target-ele-pro

        hists.computeIfAbsent("th.lt.2/hepx",hmm2epx).fill(epx.mass2())
        pi0s.each{ig1,ig2,gmass->
          def g1 = LorentzVector.withPID(22,*['px','py','pz'].collect{partb.getFloat(it,ig1)})
          def g2 = LorentzVector.withPID(22,*['px','py','pz'].collect{partb.getFloat(it,ig2)})
          def gg = g1+g2
          def misse = epx-gg
          hists.computeIfAbsent("th.lt.2/hmisse",hmm2epx).fill(misse.e())
        }
      }

      eps.collectMany{it[2]}.findAll().unique().each{ig1,ig2,ggmass->
        hists.computeIfAbsent("th.lt.2/hggm",hggm).fill(ggmass)
      }
    }
  }
}
