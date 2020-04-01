#!/home/kenjo/.groovy/coatjava/bin/run-groovy
import org.jlab.io.hipo.HipoDataSource
import org.jlab.detector.base.DetectorType
import org.jlab.clas.physics.Vector3
import org.jlab.groot.data.H1F
import org.jlab.groot.data.H2F
import org.jlab.groot.data.TDirectory
import groovyx.gpars.GParsPool
import java.util.concurrent.ConcurrentHashMap
import java.util.concurrent.Executors
import java.util.concurrent.ScheduledExecutorService
import java.util.concurrent.TimeUnit
import java.util.concurrent.atomic.AtomicInteger

MyMods.enable()
/////////////////////////////////////////////////////////////////////

def outname = args[0].split('/')[-1]

def timeline = new ConcurrentHashMap()

GParsPool.withPool 12, {
  args.eachParallel{fname->
    def reader = new HipoDataSource()
    reader.open(fname)

    //while(reader.hasEvent() && evcount.get()<5000000) {
    while(reader.hasEvent()) {
      def event = reader.getNextEvent()
      if(event.hasBank("RUN::config")) {
        def cnfb = event.getBank("RUN::config")
        def ts = cnfb.getLong("timestamp", 0)

        if(event.hasBank("REC::Particle")) {
          def partb = event.getBank("REC::Particle")
          def iele = (0..<partb.rows()).find{partb.getInt('pid',it)==11 && partb.getShort("status",it)<0}
          if(iele!=null) timeline[ts] = [1,1]
        }

        if(event.hasBank("RUN::scaler")) {
          def scaler = event.getBank("RUN::scaler")
          timeline[ts] = [scaler.getFloat('fcup',0), scaler.getFloat('fcupgated',0)]
        }
      }    
    }

    reader.close()
  }
}


def fc0 = null, fc1 = null, ts0 = null
def nelemap = [:].withDefault{0}

timeline.sort().each{ts,fcs->
  if(fcs[0]==1 && fc0!=null) {
    nelemap[[ts0,fc0,fc1]]++
  } else if(fcs[0]!=1) {
    ts0 = ts
    fc0 = fcs[0]
    fc1 = fcs[1]
  }
}

def heyield = new H1F("heyield", "ele yield", 100,0,50)
def hcurr = new H1F("hcurr", "current", 120,0,60)
def hyieldcurr = new H2F("hyieldcurr", "yeild vs current", 120,0,60, 100,0,50)

def ecounts = nelemap.sort{it.key[0]}.collect{[*it.key, it.value]}
(0..ecounts.size-2).each{i->
  def (dt,dfc0,dfc1,nele) = [*[0,1,2].collect{ecounts[i+1][it]-ecounts[i][it]}, ecounts[i][3]]

  def yld = nele/dfc1
  def curr = dfc1/dt/4*1e9

  heyield.fill(yld)
  hcurr.fill(curr)
  hyieldcurr.fill(curr, yld)
}

def out = new TDirectory()
out.mkdir("/root")
out.cd("/root")
out.addDataSet(heyield)
out.addDataSet(hcurr)
out.addDataSet(hyieldcurr)
out.writeFile("ecounts_${outname}")
