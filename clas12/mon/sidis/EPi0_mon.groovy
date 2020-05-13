package clas12.mon.sidis

import org.jlab.clas.physics.LorentzVector
import org.jlab.clas.physics.Vector3
import org.jlab.groot.data.H1F
import org.jlab.groot.data.H2F
import java.util.concurrent.ConcurrentHashMap


class EPi0_mon {
  def hists = new ConcurrentHashMap()

  def beam = LorentzVector.withPID(11,0,0,10.6)
  def target = LorentzVector.withPID(2212,0,0,0)

  def hww = {new H1F("$it","$it",200,0,5)}
  def hggm = {new H1F("$it","$it",200,0.05,0.2)}
  def hthepe = {new H2F("$it","$it",200,0,11,200,0,40)}

  def banknames = ['REC::Particle','REC::Calorimeter']
  def processEvent(event) {
    if(banknames.every{event.hasBank(it)}) {
      def (partb,ecb) = banknames.collect{event.getBank(it)}

      (0..<partb.rows()).find{partb.getInt('pid',it)==11 && partb.getShort("status",it)<0}?.with{iele->
        def chi2pid = partb.getFloat('chi2pid',iele)
        def ele = LorentzVector.withPID(11,*['px','py','pz'].collect{partb.getFloat(it,iele)})
        def ex = beam+target-ele
        def yy = (beam.e()-ele.e())/beam.e()

        def igs = (0..<partb.rows()).findAll{ig->partb.getInt('pid',ig)==22 && partb.getShort('status',ig)>=2000 &&
          ['px','py','pz'].sum{partb.getFloat(it,ig)**2}>0.36 &&
          ele.vect().theta(new Vector3(*['px','py','pz'].collect{partb.getFloat(it,ig)})) > 8
        }

        def ggs = [igs,igs].combinations().findAll{it[1]>it[0]}.each{ig1,ig2->
          def g1 = LorentzVector.withPID(22,*['px','py','pz'].collect{partb.getFloat(it,ig1)})
          def g2 = LorentzVector.withPID(22,*['px','py','pz'].collect{partb.getFloat(it,ig2)})
          def gg = g1+g2
          def ggmass = gg.mass()
          def epi0x = ex - gg
          hists.computeIfAbsent("hggm",hggm).fill(ggmass)
          if(ggmass<0.16 && ggmass>0.1 && gg.p()>1.25 && epi0x.mass()>1.5 && ex.mass()>2) {
            hists.computeIfAbsent("hggm1",hggm).fill(ggmass)
            hists.computeIfAbsent("hw",hww).fill(ex.mass())
            [1,0.85,0.8,0.75].collect{["y.lt.$it",yy<it]}.findAll{it[1]}.each{name,cut->
              hists.computeIfAbsent("$name/hthepe",hthepe).fill(ele.e(), Math.toDegrees(ele.theta()))
              if(chi2pid<3)
                hists.computeIfAbsent("chi2pid.lt.3/$name/hthepe",hthepe).fill(ele.e(), Math.toDegrees(ele.theta()))
            }
          }
        }
      }
    }
  }
}
