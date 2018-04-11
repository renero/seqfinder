import org.apache.spark.sql.SQLContext
import org.apache.spark.sql.types.TimestampType
import org.apache.spark.mllib.linalg.distributed.{CoordinateMatrix, MatrixEntry}
import org.apache.spark.HashPartitioner
import org.apache.spark.sql.SQLContext
import org.apache.spark.sql.DataFrame
import org.apache.spark.sql.functions._
import org.apache.spark.sql.functions.udf
import org.joda.time.format._
import org.joda.time._
import org.apache.spark.rdd
import org.apache.spark.rdd.RDD
import java.sql.Timestamp

val sqlContext = spark.sqlContext
import sqlContext.implicits._

case class Event(eventType: String, timestamp: Long, idOpVentas:String, desEstado:String, tipificacion:String, nps: Integer, tipoSiniestro:String, tipoReclamacion:String)
case class Touchpoint(uid: Int, event: Event)

// COMMAND ----------

def readCsv(inputFile:String, sep:String=","): DataFrame = {
  val csv = sqlContext.read.format("com.databricks.spark.csv").
      option("header", "true").
      option("delimiter", sep).
      option("dateFormat", "yyyy-MM-dd'T'HH:mm:ssZ").
      option("inferSchema", "true").
      option("nullValue", "null").
      load(inputFile)
  return csv
}

// Given a datetime string, it returns the timestamp in days. 
// We've three different formats for dates in the file...
def strToTime(s: String):Long = { 
  if(s.length == 19)
    DateTimeFormat.forPattern("yyyy-MM-dd HH:mm:ss").
      parseDateTime(s).getMillis()///(1000*60*60*24)
  else if (s.length == 8) 
     DateTimeFormat.forPattern("yyyy-MM-dd HH:mm:ss").
      parseDateTime(s.replaceAll("/","-")+(" 12:00:00")).
      getMillis()///(1000*60*60*24)
 else
    DateTime.parse(s).getMillis()///(1000*60*60*24)
}

// Given a timestamp, it returns that milliseconds figure in days since 1970
def timestampToDays(ts: Long):Long = ts/(1000*60*60*24)


// Given two timestamps it computes the difference between them both in months and returns
// a string, starting with letter 'M', and followed by the hex code for the count of months
// between 0 and 15 (max).
def encodeTimeGap(ts1:Long, ts2:Long, encodingChar:Char = 'Z'): String = {
  val numMonths = ((timestampToDays(ts2) - timestampToDays(ts1))/30)
  val gapLen = if(numMonths > 15) 15 else numMonths
  encodingChar + gapLen.toHexString
}

// Add events representing time with NO activity.
def addTimegaps(transformed: RDD[Array[Touchpoint]]): RDD[Touchpoint] = {
  val windows = transformed.filter(_.length > 1).map(userTouchPoints => userTouchPoints.sliding(2).toArray)
  val diffs = windows.map(pairs => pairs.map { 
      case Array(tp1, tp2) => Touchpoint(tp1.uid, Event(encodeTimeGap(tp1.event.timestamp, tp2.event.timestamp), tp1.event.timestamp+1, "","","",0,"",""))
    })
  val mutePeriods = diffs.map(user => user.filter(line => line.event.eventType != "Z0")).filter( !_.isEmpty ).flatMap(x => x)
  val allEvents = transformed.flatMap(x => x).union(mutePeriods)
  
  return allEvents.sortBy(tp => (tp.uid, tp.event.timestamp))
}

def mapToRDD(csv: DataFrame): RDD[Array[Touchpoint]] = 
  csv.rdd.
    map(row => Touchpoint(
      row.getAs[Int]("ID_CLIENTE"),
      Event(row.getAs[String]("ID_TOUCHPOINT"),
        strToTime(row.getAs[String]("FEC_CONTROL")),
        row.getAs[String]("ID_OPERACION_VENTAS"),
        row.getAs[String]("DES_ESTADO"),
        row.getAs[String]("DESC_TIPIFICACION_N1_GC"),
        row.getAs[Int]("NUM_NOTA_SATISFACCION"),
        row.getAs[String]("DES_TIPO_SINIESTRO"),
        row.getAs[String]("TIPO_RECLAMACION")
      )
    )).
    groupBy(tp => tp.uid).
    map(tps => tps._2.toArray.sortBy(tp => (tp.uid, tp.event.timestamp)))
    //map(tps => tps.map(tp => Touchpoint(tp.uid, Event(tp.event.eventType, timestampToDays(tp.event.timestamp), 
    //  tp.event.idOpVentas, tp.event.desEstado, tp.event.tipificacion, tp.event.nps, tp.event.tipoSiniestro, tp.event.tipoReclamacion))))


def letterEncode(touchpoint:String, operacion:String, estado:String, tipificacion:String, 
  nota: Int, siniestro: String, reclamacion: String) : String = {
  touchpoint match {
    case "PRESTACIONES" => {
      siniestro match {
        case "Reembolso Internacional" => "SI"
        case "Compensaciones Asociados Nacionales" => "Sn"
        case "Enfermedad" => "SE"
        case "Compensaciones Asociados Internacionales" => "Si"
        case "Reembolso Nacional" => "SN"
        case _ => "__"
      }
    }
    case "RECLAMACIONES" => {
      reclamacion match {
        case "EUROP" => "Re"
        case "SERVICIOS GENERALES" => "Rg"
        case "PROMOCION" => "Rp"
        case "LEY PROTECCION DATOS" => "Rl"
        case "CONSULTA LOPD" => "Rq"
        case "WELCOME PACK" => "RW"
        case "SUBIDA DE PRIMAS" => "RS"
        case "COPAGO" => "RC"
        case "AGRADECIMIENTO" => "Rt"
        case "SERVICIOS FINANCIEROS" => "RF"
        case "SUGERENCIAS" => "Rs"
        case "ATENCION ASISTENCIAL" => "RA"
        case "ATENCION ADMINISTRATIVA" => "Ra"
        case "BAJAS" => "RB"
        case "GESTION DE POLIZAS" => "RG"
        case "CONDICIONADO GENERAL/COBERTURA" => "Rc"
        case "TARJETA" => "RT"
        case "GUIA MEDICA" => "RM"
        case "WEB" => "Rw"
        case "VOLANTES" => "Rv"
        case "INFORMACION" => "Ri"
        case _ => "__"
      }    
    }
    case "AUTORIZACIONES" => {
      estado match {
        case "AUTORIZADO" => "AA"
        case "DENEGADO" =>   "AD"
        case "PENDIENTE" =>  "AP"
        case _ => "__"
      }
    }
    case "MOVIMIENTOS" => {
      operacion match {
          case "B" =>  "BB"
          case "TI" => "MA"
          case "TO" => "MM"
          case "A" =>  "AL"
          case _ => "__"
      }
    }
    case "RECIBOS" => {
      estado match {
        case "ABIERTO" => "FO"
        case "ANULADO" => "Fc"
        case "ANULADO POR BAJA DE POLIZA" => "FR"
        case "AUTORIZADO" => "FA"
        case "CERRADO" => "Fx"
        case "COBRADO" => "FC"
        case "COBRO ANULADO" => "FN"
        case "COMPENSADO" => "FB"
        case "DENEGADO" => "FD"
        case "IMPAGADO" => "FI"
        case "PAGADO" => "FP"
        case "PENDIENTE" => "Fp"
        case _ => "__"
      }
    }
    case str if str.startsWith("NPS A") => {
      if(nota <= 5) "NN" else "NP"
    }
    case str if str.startsWith("NPS N") => {
      if(nota <= 5) "Nn" else "Np"
    }
    case "GESTOR DE CONTACTOS" => {
      tipificacion match {
        case s if s.startsWith("ACTUALIZACI") => "GU"
        case s if s.startsWith("AUTORIZACI") => "GA"
        case s if s.startsWith("BAJA") => "GB"
        case s if s.startsWith("CAMPA") => "Gt"
        case s if s.startsWith("CITA") => "Gp"
        case s if s.startsWith("CONSULTA BENEFICIARIO") => "GQ"
        case s if s.startsWith("CONTRATACI") => "GH"
        case s if s.startsWith("COPAGO") => "Gc"
        case s if s.startsWith("DENTAL") => "GD"
        case s if s.startsWith("DOCUMENTACI") => "Gd"
        case s if s.startsWith("EMISI") => "Ge"
        case s if s.startsWith("FACTURACI") => "GF"
        case s if s.startsWith("FIDELIZACI") => "Gf"
        case s if s.startsWith("GRACIABLE") => "Gg"
        case s if s.startsWith("INCIDENCIA") => "GI"
        case s if s.startsWith("INFORMACI") => "Gi"
        case s if s.startsWith("LLAMADA") => "GL"
        case s if s.startsWith("MI SANITAS/APP") => "Ga"
        case s if s.startsWith("MODIFICACI") => "Gu"
        case s if s.startsWith("RECLAMACI") => "GR"
        case s if s.startsWith("REEMBOLSO") => "Gr"
        case s if s.startsWith("SEGUIMIENTO COMUNICACIONES MKT") => "Gm"
        case s if s.startsWith("TIS") => "GT"
        case s if s.startsWith("WEB") => "Gw"
        case _ => "__"
      }
    }
    case "DIGITAL" => "DD"
    case "ALTA" => "AA"
    case "TP FICITICIO" => "__"
    case "GASTO MEDICO" => "GM"
    case "COPAGOS" => "CP"
    case "OOP DENTAL" => "OD"
    case "BAJA" => "BB"
    case "ACTUALIZACION" => "MA"
    case "MIGRACION" => "MM"
    case "Z1" => "Z1"
    case "Z2" => "Z2"
    case "Z3" => "Z3"
    case "Z4" => "Z4"
    case "Z5" => "Z5"
    case "Z6" => "Z6"
    case "Z7" => "Z7"
    case "Z8" => "Z8"
    case "Z9" => "Z9"
    case "Za" => "Za"
    case "Zb" => "Zb"
    case "Zc" => "Zc"
    case "Zd" => "Zd"
    case "Ze" => "Ze"
    case "Zf" => "Zf"
    case _ => "__"
  }
}


// COMMAND ----------

def encode(data: DataFrame): RDD[(Int, String, Long, String)] = {
  // Function to re-encode the MOVIMIENTO TO correctly reflect "churn".
  // Reduce dataset
  val reduced = data.
    select($"ID_CLIENTE",$"ID_TOUCHPOINT",$"FEC_CONTROL", $"ID_OPERACION_VENTAS", $"DES_ESTADO", 
      $"DESC_TIPIFICACION_N1_GC", $"NUM_NOTA_SATISFACCION", $"DES_TIPO_SINIESTRO", $"TIPO_RECLAMACION").
      filter("FEC_CONTROL is not null")

  // Inject timegaps
  val enriched = addTimegaps(mapToRDD(reduced)).
      map {
        case touchpoint =>
          (touchpoint.uid, touchpoint.event.eventType, touchpoint.event.timestamp, 
            letterEncode(touchpoint.event.eventType, touchpoint.event.idOpVentas, touchpoint.event.desEstado, touchpoint.event.tipificacion, 
              touchpoint.event.nps, touchpoint.event.tipoSiniestro, touchpoint.event.tipoReclamacion))
      }

    return enriched.filter(element => element._1 != 0)
}

// COMMAND ----------

def buildSeq(rdd: RDD[(Int, String, Long, String)], remove_repetitions:Boolean = true): RDD[(String)] = {
  val seqs = rdd.map {
      case (uid, touchpoint, timstamp, event) => (uid, event)
    }.
    groupByKey.
    map(tuple => (tuple._1, tuple._2.mkString("")))

  // Check if we need to take repetitions off.
  // "(.)\\1{1,}" means any character (added to group 1) followed by itself at least once
  // "$1" references contents of group 1
  if(remove_repetitions)
    seqs.map(seq => (seq._1, seq._2.replaceAll("(..)\\1{1,}", "$1"))).  
      map(seq => seq._1+","+seq._2)
  else
    seqs.map(seq => seq._1+","+seq._2)
}

// COMMAND ----------

def process(csvFile:String, outputFile:String, sep: String = ",") {
    val csv  = readCsv(csvFile, sep)
    val data = encode(csv)
    val seqs = buildSeq(data)  
    seqs.coalesce(1).saveAsTextFile(outputFile)
}

// COMMAND ----------
val dataPath = "adl://dlakepoctestadls.azuredatalakestore.net/sanitas/incoming/tablon_bajas/"
val altas = "adl://dlakepoctestadls.azuredatalakestore.net/sanitas/incoming/tablon_bajas/tablon_altas_v3_unquoted_ncnf/tablon_altas_v3_unquoted_ncnf.csv.bz2"
val bajas = "adl://dlakepoctestadls.azuredatalakestore.net/sanitas/incoming/tablon_bajas/tablon_bajas_v3_unquoted_ncnf/tablon_bajas_v3_unquoted_ncnf.csv.bz2"

val bajasFile  = readCsv(bajas, ",")
val dataB = encode(bajasFile)
val seqsB = buildSeq(dataB)  
seqsB.coalesce(1).saveAsTextFile(dataPath + "seqs_bajas_ncnf_norep_v4.0")

val altasFile = readCsv(altas, ",")
val dataA = encode(altasFile)
val seqsA = buildSeq(dataA)  
seqsA.coalesce(1).saveAsTextFile(dataPath + "seqs_altas_ncnf_norep_v4.0")

