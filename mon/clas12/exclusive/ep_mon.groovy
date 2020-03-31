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


class ep_mon {
  def hists = new ConcurrentHashMap()

  def beam = LorentzVector.withPID(11,0,0,10.6)
  def target = LorentzVector.withPID(2212,0,0,0)

  def hm2pro = {new H1F("$it","$it",200,0,2)}
  def hdfi = {new H1F("$it","$it",200,150, 210)}

  def he0 = {new H1F("$it","$it",200,7,13)}
  def he0th = {new H2F("$it","$it",200,7,13,200,7,13)}

  def banknames = ['REC::Particle','REC::Calorimeter']
  def processEvent(event) {
    if(banknames.every{event.hasBank(it)}) {
      def (partb,ecb) = banknames.collect{event.getBank(it)}

      (0..<partb.rows()).findAll{partb.getInt('pid',it)==11 && partb.getShort("status",it)<0}
        .collectMany{iele->
          (0..<partb.rows()).findAll{partb.getInt('pid',it)==2212 && partb.getShort('status',it)>=4000}.collect{ipro->[iele,ipro]}
        }.each{iele,ipro->
          def ele = LorentzVector.withPID(11,*['px','py','pz'].collect{partb.getFloat(it,iele)})
          def pro = LorentzVector.withPID(2212,*['px','py','pz'].collect{partb.getFloat(it,ipro)})

          def esec = (0..<ecb.rows()).find{ecb.getShort('pindex',it)==0}.with{ecb.getByte('sector',it)}

          def ex = beam+target-ele
          def dfi = Math.toDegrees((ele.phi()-pro.phi()).abs())
          def e0 = 0.938*(-1+1/(Math.tan(ele.theta()/2)*Math.tan(pro.theta())))

          hists.computeIfAbsent("hW:$esec",hm2pro).fill(ex.mass2())
          hists.computeIfAbsent("hdfi:$esec",hdfi).fill(dfi)

          [
            ['dfi.lt.4', (dfi-180).abs()<4],
          ].findAll{it[1]}.each{name,cut->
            hists.computeIfAbsent("$name/hW:$esec",hm2pro).fill(ex.mass2())
            hists.computeIfAbsent("$name/hdfi:$esec",hdfi).fill(dfi)
            hists.computeIfAbsent("$name/he0:$esec",he0).fill(e0)
            hists.computeIfAbsent("$name/he0th:$esec",he0th).fill(Math.toDegrees(ele.theta()), e0)
          }
        }
    }
  }
}
