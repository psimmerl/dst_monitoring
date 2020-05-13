package clas12.mon.sidis

import org.jlab.clas.physics.LorentzVector
import org.jlab.clas.physics.Vector3
import org.jlab.groot.data.H1F
import org.jlab.groot.data.H2F
import java.util.concurrent.ConcurrentHashMap


class EPip_mon {
  def hists = new ConcurrentHashMap()

  def beam = LorentzVector.withPID(11,0,0,10.6)
  def target = LorentzVector.withPID(2212,0,0,0)

  def hww = {new H1F("$it","$it",200,0,5)}
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

        def pips = (0..<partb.rows()).findAll{ipip->partb.getInt('pid',ipip)==211 && ['px','py','pz'].sum{partb.getFloat(it,ipip)**2}>0.36}
          .findResults{ipip->
            def pip = LorentzVector.withPID(211,*['px','py','pz'].collect{partb.getFloat(it,ipip)})
            def epi0x = ex - pip
            if(ex.mass()>2 && epi0x.mass()>1.5 && pip.p()>1.25) {
              return (partb.getShort('status',ipip)/2000).toInteger()
            } else {
              return null
            }
          }

        [1,2].findAll{detid->pips.any{detid==it}}.each{detid->
           def det = detid==1 ? 'FD' : 'CD'
           hists.computeIfAbsent("hw:${det}",hww).fill(ex.mass())
           [1,0.85,0.8,0.75].collect{["y.lt.$it",yy<it]}.findAll{it[1]}.each{name,cut->
              hists.computeIfAbsent("$name/hthepe:$det",hthepe).fill(ele.e(), Math.toDegrees(ele.theta()))
              if(chi2pid<3)
                hists.computeIfAbsent("chi2pid.lt.3/$name/hthepe:$det",hthepe).fill(ele.e(), Math.toDegrees(ele.theta()))
           }
        }
      }
    }
  }
}
