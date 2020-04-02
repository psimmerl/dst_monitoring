#!/home/kenjo/.groovy/coatjava/bin/run-groovy
import org.jlab.io.hipo.HipoDataSource
import org.jlab.groot.data.H1F
import org.jlab.groot.data.H2F
import org.jlab.groot.data.TDirectory
import groovyx.gpars.GParsPool
import java.util.concurrent.ConcurrentHashMap
import java.util.concurrent.Executors
import java.util.concurrent.ScheduledExecutorService
import java.util.concurrent.TimeUnit
import java.util.concurrent.atomic.AtomicInteger
import clas12.mon.exclusive.EP_mon
import clas12.mon.exclusive.EPPi0_mon
import clas12.mon.exclusive.EPPipPim_mon
import clas12.mon.exclusive.ENPip_mon
import clas12.mon.fcup.FCup_mon
import clas12.mon.electron.ECounts
import clas12.mon.ft.FT_mon
import clas12.groovy.Sugar

Sugar.enable()
/////////////////////////////////////////////////////////////////////

def outname = args[0].split('/')[-1]

def processors = [new EPPi0_mon(), new EPPipPim_mon(), new EP_mon()]

def evcount = new AtomicInteger()
def save = {
  processors.each{
    def out = new TDirectory()
    out.mkdir("/root")
    out.cd("/root")
    it.hists.each{out.writeDataSet(it.value)}
    def clasname = it.getClass().getSimpleName()
    out.writeFile("mon_${clasname}_${outname}")
  }
  println "event count: "+evcount.get()
  evcount.set(0)
}

def exe = Executors.newScheduledThreadPool(1)
exe.scheduleWithFixedDelay(save, 5, 30, TimeUnit.SECONDS)

GParsPool.withPool 12, {
  args.eachParallel{fname->
    def reader = new HipoDataSource()
    reader.open(fname)

    while(reader.hasEvent()) {
      evcount.getAndIncrement()
      def event = reader.getNextEvent()
      processors.each{it.processEvent(event)}
    }

    reader.close()
  }
}

processors.each{if(it.metaClass.respondsTo(it, 'finish')) it.finish()}

exe.shutdown()
save()
