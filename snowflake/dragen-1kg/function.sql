-- Copyright (c) 2022 Snowflake Inc.  All Rights Reserved


-- UDTF to ingest gzipped vcf file.  If num_lines > 0, only first N data lines are processed
create or replace function ingest_vcf(vcfFile string, num_lines integer)
returns table (chrom string, pos integer, id string, ref string, alt string, qual string,
                   filter string, info variant, sampleVals variant, sampleId string)
language java
handler='App'
target_path='@~/ingest_vcf.jar'
as
$$
import java.io.*;
import java.util.*;
import java.util.zip.GZIPInputStream;
import java.util.stream.Stream;


class OutputRow {
   public String chrom; 
   public int    pos;      
   public String id; 
   public String ref; 
   public String alt;  
   public String qual;  
   public String filter;
   public String info;
   public String sampleVals;
   public String sampleId;


   public OutputRow(String chrom, int pos, String id, String ref, String alt, String qual,
                     String filter, String info, String sampleVals, String sampleId) {


       this.chrom      = chrom; 
       this.pos        = pos;   
       this.id         = id;  
       this.ref        = ref;   
       this.alt        = alt;  
       this.qual       = qual;  
       this.filter     = filter;
       this.info       = info;
       this.sampleVals = sampleVals;
       this.sampleId   = sampleId;


   }
}


class DictEntry {
   public String number;
   public String type;
   public String description; 
  
   public DictEntry(String number, String type, String description) {
       this.number = number;
       this.type = type;
       this.description = description;
   }
}


public class App {


   public static Class getOutputClass() {
       return OutputRow.class;
   }


   // Return quoted values unless numeric, to populate JSON values
   // Also replace "dots" with empty strings
   public static String delimitedValue(String rawValue, String type) {
       String theValue = (rawValue.equals(".")) ? "" : rawValue;
       return ((type.equals("Integer") || type.equals("Float")) ? theValue : "\"" + theValue + "\"");
   }


   // Assemble a Json Output object element from a VCF element
   public static String createJsonElement(String theKey, String unparsedValue, Map<String, DictEntry> dict) {


       // Lookup metadata in the dictionary
       String theNumber = dict.get(theKey).number;
       String theType = dict.get(theKey).type;


       String theValue = (theNumber.equals("0") ? "" : unparsedValue);


       if (theNumber.equals("0")) { // Boolean
           return("\"" + theKey + "\": true");
       }
       else if (theNumber.equals("1")) { // Singleton
           return("\"" + theKey + "\": " + delimitedValue(theValue, theType));
       }
       else { // Array of values
           List<String> vals = new ArrayList<String>();
               for (String e : theValue.split(",")) {
                   vals.add(delimitedValue(e, theType));
               }
           return("\"" + theKey + "\": [" + String.join(",", vals) + "]");
       }
   }


   // Parse content of Info fields into JSON
   public static String createInfoJson(Map<String, DictEntry> dict, String infoText) {
       if (infoText.equals(".")) {
           return null;
       }
       List<String> elements  = new ArrayList<String>();
       String[] infoArray = infoText.split(";");
       for (String item : infoArray) {
          String theKey = item.split("=")[0].trim();
          String unparsedValue = (item.contains("=") ? item.split("=")[1].trim() : "");
          elements.add(createJsonElement(theKey, unparsedValue, dict));
       }
       return ("{" + String.join(",", elements) + "}");
   }


   // Parse genomic values based on Format field into Json
   public static String createFormatJson(Map<String, DictEntry> dict, String formatString, String formatVals) {
       List<String> elements = new ArrayList<String>();
       String[] formatArray = formatString.split(":");  // available fields
       String[] valsArray = formatVals.split(":"); // the values -- may be truncated relative to format entries
       for (int i = 0 ; i < valsArray.length; i++ ) {
           String theKey = formatArray[i];
           String unparsedValue = valsArray[i];
           // add an entry to the JSON so long as the field is not empty (!== ".")
           if (! unparsedValue.equals(".")) {
               elements.add(createJsonElement(theKey, unparsedValue, dict));                   
           }  
       }                  
       return ("{" + String.join(",", elements) + "}");
   }


   // Parse a line from header INFO or FORMAT entries and populate the respective dictionaries
   public static void addDictionaryElement(Map<String, DictEntry> dict, String descriptor) {
       String infoData = descriptor.split("=<")[1];
       dict.put(infoData.split("ID=")[1].split(",")[0],
               new DictEntry(
                   infoData.split("Number=")[1].split(",")[0],
                   infoData.split("Type=")[1].split(",")[0],
                   infoData.split("Description=")[1].split(">$")[0]
               ));       
   }


   public static Stream<OutputRow> process(InputStream vcfFile, long num_lines) throws IOException {


       Map<String, DictEntry> infos = new HashMap<>();
       Map<String, DictEntry> formats = new HashMap<>();


       try {
           GZIPInputStream gzStream = new GZIPInputStream(vcfFile);
           Reader vcf = new InputStreamReader(gzStream, "US-ASCII");
           BufferedReader vcfBuf = new BufferedReader(vcf);
           String thisLine = null;
           while ((thisLine = vcfBuf.readLine()) != null) {
               if (! thisLine.startsWith("##")) {
                   break;
               }
               else if (thisLine.startsWith("##INFO")) {
                   addDictionaryElement(infos, thisLine);
               }
               else if (thisLine.startsWith("##FORMAT")) {
                   addDictionaryElement(formats, thisLine);
               }               
           }
          
           if (! thisLine.startsWith("#CHROM")) {
               throw new IOException("Unknown VCF Layout, expected '#CHROM' header line");
           }


           final String[] columns = thisLine.split("#")[1].split("\t");


           // Map the position of all "Well Known Columns"
           Map<String,Integer> columnIndex = new HashMap<>(10);
           String knownCols = "CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT";
           for (int i = 0; i < columns.length; i++) {
               columnIndex.put(columns[i],i);
               if (knownCols.indexOf(columns[i]) == -1) {
                   break;
               }
           }


           // 0-based index of first sample positon for single or multisample VCF
           // Add some validation...
           final int firstSampleCol = (columns.length > 9) ? 9 :
           (columns.length == 9 && columnIndex.get("FORMAT") == null) ? 8 : -1;


           // Handle the data lines
          
           return ((num_lines == 0) ? vcfBuf.lines() : vcfBuf.lines().limit(num_lines)).flatMap( dataLine -> {
               String[] vals = dataLine.split("\t");
               // Get the sample-independent data
               String chrom  =  vals[columnIndex.get("CHROM")];
               int pos       =  Integer.parseInt(vals[columnIndex.get("POS")]);
               String id     =  vals[columnIndex.get("ID")];
               String ref    =  vals[columnIndex.get("REF")];
               String alt    =  vals[columnIndex.get("ALT")];     
               String qual   =  vals[columnIndex.get("QUAL")];  
               String filter =  vals[columnIndex.get("FILTER")];  
               String info   =  createInfoJson(infos, vals[columnIndex.get("INFO")]);
               return java.util.stream.IntStream.range((firstSampleCol >= 0 ? firstSampleCol : columns.length - 1), columns.length)
                   .mapToObj( i -> {
                       String format =  (columnIndex.get("FORMAT") == null) ? null : vals[columnIndex.get("FORMAT")];
                       String sampleVals = (format == null) ?      // case where there is no FORMAT column,
                                               null :              // e.g. no populated per-sample column of data in the VCF
                                           (firstSampleCol >= 0) ? // case where there is at least one populated sample column
                                           createFormatJson(formats, format, vals[i]) :
                                               null;               // case where no sample column is defined at all (annotation VCF)
                       String sampleId   = (firstSampleCol >= 0) ? columns[i] : "";


                       return new OutputRow(chrom, pos, id, ref, alt, qual, filter, info, sampleVals, sampleId);     
                   });
           });
       }
       finally{}
   }
}
$$
;


-- UDFs to help parse out allele values and reduce query complexity!
CREATE OR REPLACE SECURE FUNCTION allele1 (ref varchar, alt array, gt varchar)
RETURNS varchar as
$$
case when contains(gt, '|') then
   case when strtok(gt, '|', 1) = '0' then ref else alt[strtok(gt, '|', 1)::integer-1] end
else
   case when strtok(gt, '/', 1) = '0' then ref else alt[strtok(gt, '/', 1)::integer-1] end
end
$$
;


CREATE OR REPLACE SECURE FUNCTION allele2 (ref varchar, alt array, gt varchar)
RETURNS varchar as
$$
case when contains(gt, '|') then
   case when strtok(gt, '|', 2) = '0' then ref else alt[strtok(gt, '|', 2)::integer-1] end
else
   case when strtok(gt, '/', 2) = '0' then ref else alt[strtok(gt, '/', 2)::integer-1] end
end
$$
;