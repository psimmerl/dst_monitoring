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


class epippim_mon {
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
          (0..<partb.rows()).findAll{partb.getInt('pid',it)==211}.collect{ipip->[iele,ipip]}
        }.collectMany{iele,ipip->
          (0..<partb.rows()).findAll{partb.getInt('pid',it)==-211}.collect{ipim->[iele,ipip,ipim]}
        }.each{iele,ipip,ipim->
          def ele = LorentzVector.withPID(11,*['px','py','pz'].collect{partb.getFloat(it,iele)})
          def pip = LorentzVector.withPID(211,*['px','py','pz'].collect{partb.getFloat(it,ipip)})
          def pim = LorentzVector.withPID(-211,*['px','py','pz'].collect{partb.getFloat(it,ipim)})
          def impippim = pip+pim

          hists.computeIfAbsent("impippim",him).fill(impippim.mass())
        }
    }
  }
}
