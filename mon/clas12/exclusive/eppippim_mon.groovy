package mon.clas12.exclusive

import org.jlab.detector.base.DetectorType
import org.jlab.clas.physics.LorentzVector
import org.jlab.clas.physics.Vector3
import org.jlab.groot.data.H1F
import org.jlab.groot.data.H2F
import java.util.concurrent.ConcurrentHashMap
import java.util.concurrent.ConcurrentLinkedQueue
import java.util.concurrent.atomic.AtomicInteger
import java.util.concurrent.Executors
import java.util.concurrent.ScheduledExecutorService
import java.util.concurrent.TimeUnit
import java.nio.ByteBuffer


class eppippim_mon {
  def hists = new ConcurrentHashMap()

  def beam = LorentzVector.withPID(11,0,0,10.6)
  def target = LorentzVector.withPID(2212,0,0,0)

  def hmm2 = {new H1F("$it","$it",250,-0.5,2)}
  def him = {new H1F("$it","$it",180,0.2,2)}

  def banknames = ['REC::Particle']
  def processEvent(event) {
    if(banknames.every{event.hasBank(it)}) {
      def (partb) = banknames.collect{event.getBank(it)}

      (0..<partb.rows()).findAll{partb.getInt('pid',it)==11 && partb.getShort("status",it)<0}
        .collectMany{iele->
          (0..<partb.rows()).findAll{partb.getInt('pid',it)==2212}.collect{ipro->[iele,ipro]}
        }.collectMany{iele,ipro->
          (0..<partb.rows()).findAll{partb.getInt('pid',it)==211}.collect{ipip->[iele,ipro,ipip]}
        }.collectMany{iele,ipro,ipip->
          (0..<partb.rows()).findAll{partb.getInt('pid',it)==-211}.collect{ipim->[iele,ipro,ipip,ipim]}
        }.each{iele,ipro,ipip,ipim->
          def ele = LorentzVector.withPID(11,*['px','py','pz'].collect{partb.getFloat(it,iele)})
          def pro = LorentzVector.withPID(2212,*['px','py','pz'].collect{partb.getFloat(it,ipro)})
          def pip = LorentzVector.withPID(211,*['px','py','pz'].collect{partb.getFloat(it,ipip)})
          def pim = LorentzVector.withPID(-211,*['px','py','pz'].collect{partb.getFloat(it,ipim)})
          def mm2pip = beam+target-ele-pro-pim
          def mm2pim = beam+target-ele-pro-pip
          def mm2pro = beam+target-ele-pip-pim
          def impropip = pro+pip
          def impropim = pro+pim
          def impippim = pip+pim

          def prodet = (partb.getShort('status',ipro)/1000).toInteger()==2 ? 'FD':'CD'
          def pipdet = (partb.getShort('status',ipip)/1000).toInteger()==2 ? 'FD':'CD'
          def pimdet = (partb.getShort('status',ipim)/1000).toInteger()==2 ? 'FD':'CD'

          hists.computeIfAbsent("mm2pip:$pimdet",hmm2).fill(mm2pip.mass2())
          hists.computeIfAbsent("mm2pim:$pipdet",hmm2).fill(mm2pim.mass2())
          hists.computeIfAbsent("mm2pro:$prodet",hmm2).fill(mm2pro.mass2())
          hists.computeIfAbsent("impropip",him).fill(impropip.mass())
          hists.computeIfAbsent("impropim",him).fill(impropim.mass())
          hists.computeIfAbsent("impippim",him).fill(impippim.mass())

          [
            ['m2pip.lt.0.5/m2pim.lt.0.5', mm2pim.mass2().abs()<0.5 && mm2pip.mass2().abs()<0.5],
          ].findAll{it[1]}.each{name,cut->
            hists.computeIfAbsent("$name/mm2pip:$pimdet",hmm2).fill(mm2pip.mass2())
            hists.computeIfAbsent("$name/mm2pim:$pipdet",hmm2).fill(mm2pim.mass2())
            hists.computeIfAbsent("$name/mm2pro:$prodet",hmm2).fill(mm2pro.mass2())
            hists.computeIfAbsent("$name/impropip",him).fill(impropip.mass())
            hists.computeIfAbsent("$name/impropim",him).fill(impropim.mass())
            hists.computeIfAbsent("$name/impippim",him).fill(impippim.mass())
          }
        }
    }
  }
}
