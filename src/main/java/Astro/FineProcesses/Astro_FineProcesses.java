package Astro.FineProcesses;


import ij.IJ;
import ij.ImagePlus;
import ij.gui.Roi;
import ij.Prefs;
import ij.gui.Plot;
import ij.gui.WaitForUserDialog;
import ij.io.FileSaver;
import ij.measure.Calibration;
import ij.measure.ResultsTable;
import ij.plugin.Duplicator;
import ij.plugin.PlugIn;
import ij.plugin.filter.RankFilters;
import ij.plugin.frame.RoiManager;
import ij.process.ImageProcessor;
import ij.process.ImageStatistics;
import ij.plugin.ImageCalculator;
import java.awt.Color;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.logging.Level;
import java.util.logging.Logger;
import loci.common.services.DependencyException;
import loci.common.services.ServiceException;
import loci.common.services.ServiceFactory;
import loci.formats.FormatException;
import loci.formats.meta.IMetadata;
import loci.formats.services.OMEXMLService;
import loci.plugins.in.ImporterOptions;
import loci.plugins.util.ImageProcessorReader;
import mcib3d.geom.Point3D;
import mcib3d.image3d.ImageHandler;
import org.apache.commons.io.FilenameUtils;
import sc.fiji.localThickness.Clean_Up_Local_Thickness;
import sc.fiji.localThickness.EDT_S1D;
import sc.fiji.localThickness.Local_Thickness_Parallel;
import sc.fiji.analyzeSkeleton.AnalyzeSkeleton_;
import sc.fiji.analyzeSkeleton.SkeletonResult;
import sc.fiji.analyzeSkeleton.Graph;
import sc.fiji.analyzeSkeleton.Vertex;
import sc.fiji.analyzeSkeleton.Point;
import mcib3d.utils.ArrayUtil;



/**
 *
 * @author phm
 */



public class Astro_FineProcesses implements PlugIn {
    
    // Image directory
    public static String imageDir;
    private final boolean canceled = false;
    // output directory
    public static String outDirResults = "";
    public static Calibration cal = new Calibration();
    
    
    /*Median filter 
     * 
     * @param img
     * @param size
     */ 
    public static void median_filter(ImagePlus img, double size) {
        RankFilters median = new RankFilters();
        for (int s = 1; s <= img.getNSlices(); s++) {
            img.setZ(s);
            median.rank(img.getProcessor(), size, RankFilters.MEDIAN);
            img.updateAndDraw();
        }
    }
    
    /**
     * Flush and close images
     * @param img 
     */
    public static void flush_close(ImagePlus img) {
        img.flush();
        img.close();
    } 
    
     /**
     * Find Z with max intensity in stack
     * @param img
     * @return z
     */
    
    private static int find_max(ImagePlus img) {
        double max = 0;
        int zmax = 0;
        for (int z = 1; z <= img.getNSlices(); z++) {
            ImageProcessor ip = img.getStack().getProcessor(z);
            ImageStatistics statistics = new ImageStatistics().getStatistics(ip, ImageStatistics.MEAN, img.getCalibration());
            double meanInt = statistics.mean;
            if (meanInt > max) {
                max = meanInt;
                zmax = z;
            }
        }
        return(zmax);
    }
    
    /**
     * Threshold images, fill holes
     * @param img
     * @param thMed
     * @param fill 
     * @param calculate 
     */
    public static void threshold(ImagePlus img, String thMed, boolean fill, boolean calculate) {
        //  Threshold and binarize
        img.setZ(find_max(img));
        img.updateAndDraw();
        IJ.setAutoThreshold(img, thMed + " dark");
        Prefs.blackBackground = false;
        String method = "";
        if (calculate)
            method = "method="+thMed+" background=Dark calculate";
        else
            method = "method="+thMed+" background=Dark";
        IJ.run(img, "Convert to Mask", method);
        if (fill) {
            IJ.run(img,"Fill Holes", "stack");
        }
    }
    
    /**
    * Find background image intensity
    * find mean intensity outside astrocyte mask in astrocyte image
    * don't take in account pixels with 0 value
    * return mean intensity
     * @param imgAstro
    */
    public static double find_background(ImagePlus imgAstro) {
        ImagePlus img = new Duplicator().run(imgAstro, 3,imgAstro.getNSlices());
        median_filter(img, 0.5);
        ImagePlus imgMask = img.duplicate();
        threshold(imgMask, "Li", false, true);
        IJ.run(imgMask, "Invert", "stack");
        ImageCalculator imgCal = new ImageCalculator();
        ImagePlus img1 = imgCal.run("Multiply create stack", img, imgMask);
        IJ.run(img1, "Divide...", "value=255 stack");
        double sectionInt = 0;
        for (int n = 1; n <= img1.getNSlices(); n++) {
            double pixelInt = 0;
            int voxelZero = 0;
            img1.setZ(n);
            ImageProcessor ip = img1.getProcessor();
            for (int x = 0; x < img1.getHeight(); x++) {
                for (int y = 0; y < img1.getWidth(); y++) {
                    double voxelInt = ip.getPixelValue(x, y);
                    if ( voxelInt > 0) {
                        pixelInt += voxelInt;
                    }
                    else {
                        voxelZero++;
                    }
                }
            }
            sectionInt += pixelInt/((img1.getHeight()*img1.getWidth()) - voxelZero);
        } 
        double bg = sectionInt/img1.getNSlices();
        System.out.println("Background = "+bg);
        flush_close(img);
        flush_close(imgMask);
        flush_close(img1);
        return(bg);
    }
    /**
     * compute local thickness
     * @param imgAstro
     * @return astroMap
    **/
    public static ImagePlus localThickness3D (ImagePlus imgAstro) {
        EDT_S1D edt = new EDT_S1D();
        edt.runSilent = true;
        edt.thresh = 1;
        edt.inverse = false;
        edt.showOptions = false;
        edt.setup("", imgAstro);
        edt.run(imgAstro.getProcessor());
        ImagePlus imgEDT = edt.getResultImage();
        imgEDT.setCalibration(cal);
        Local_Thickness_Parallel locThk = new Local_Thickness_Parallel();
        locThk.runSilent = true;
        locThk.setup("", imgEDT);
        locThk.run(imgEDT.getProcessor());
        ImagePlus imgLocThk = locThk.getResultImage();
        imgLocThk.setCalibration(cal);
        Clean_Up_Local_Thickness cleanUp = new Clean_Up_Local_Thickness();
        cleanUp.runSilent = true;
        cleanUp.setup("", imgLocThk);
        cleanUp.run(imgLocThk.getProcessor());
        ImagePlus astroMap = cleanUp.getResultImage();
        // calibrate intensity to µm
        cal.setValueUnit("µm");
        astroMap.setCalibration(cal);
        IJ.run(astroMap, "Multiply...", "value="+cal.pixelWidth+" stack");
        flush_close(imgEDT);
        flush_close(imgLocThk);
        return(astroMap);
    }
    
    
    private ImagePlus findMask(ImagePlus img, Roi roi, double bg) {
        ImagePlus imgAstroMask = img.duplicate();
        threshold(imgAstroMask, "Triangle", false, false);
        imgAstroMask.setRoi(roi);
        IJ.run("Colors...", "foreground=black background=white selection=yellow");
        IJ.run(imgAstroMask, "Clear Outside","stack");
        imgAstroMask.deleteRoi();
        return(imgAstroMask);
    }
    
    
     // For each skeleton calculate diameter using image map
    public static ResultsTable diameterAnalysis (ImagePlus imgSkel, ImagePlus imgMap, String rootName, int r) throws IOException {
        ImageHandler imh = ImageHandler.wrap(imgMap);
        AnalyzeSkeleton_ analyzeSkeleton = new AnalyzeSkeleton_();
        AnalyzeSkeleton_.calculateShortestPath = true;
        analyzeSkeleton.setup("",imgSkel);
        SkeletonResult skeletonResults = analyzeSkeleton.run(AnalyzeSkeleton_.NONE, false, true, null, true, false);
        int[] branchNumbers = skeletonResults.getBranches();
        double[] branchLengths = skeletonResults.getAverageBranchLength();
        int[] nbEndPoints = skeletonResults.getEndPoints();
        int[] junctions = skeletonResults.getJunctions();
        // find points in skeleton
        Graph[] graphs = skeletonResults.getGraph();
        ArrayList<Double> diameter = new ArrayList();
        for (Graph graph : graphs) {
            ArrayList<Vertex> vertices = graph.getVertices();
            for (Vertex vertice : vertices) {
                ArrayList<Point> pts = vertice.getPoints();
                for (Point pt : pts) {
                    Point3D pt3D = new Point3D(pt.x, pt.y, pt.z);
                    diameter.add((double)imh.getPixel(pt3D));
                }
            }
        }
        // Plot diameters
        ArrayUtil diameters = new ArrayUtil(diameter);
        double mean = diameters.getMean();
        double sd = diameters.getStdDev();
        DecimalFormat dec = new DecimalFormat("#0.0000");
        ArrayUtil[] histo = diameters.getHistogram(100);
        Plot dPlot = new Plot("Astrocyte diameters plot", "", "");
        dPlot.setLimits(0, histo[0].getMaximum(), 0, histo[1].getMaximum()+10);
        dPlot.setColor(Color.black, Color.black);
        dPlot.add("bar", histo[0].getArray(), histo[1].getArray());
        dPlot.addLabel(0.4, 0.2, "Mean diameter = "+ dec.format(mean) + " +- " + dec.format(sd));
        dPlot.draw();
        dPlot.addToStack();
        
        // Plot branch Numbers
        //Plot dPlotSkel = new Plot("Astrocyte branchs plot", "Branchs", "");
        ArrayUtil branchs = new ArrayUtil(branchNumbers);
        double meanBranchs = branchs.getMean();
        double sdBranchs = branchs.getStdDev();
        ArrayUtil[] histoBranchs = branchs.getHistogram(100);
        dPlot.setColor(Color.blue, Color.blue);
        dPlot.add("bar", histoBranchs[0].getArray(), histoBranchs[1].getArray());
        dPlot.addLabel(0.4, 0.2, "Mean branchs number = "+ dec.format(meanBranchs) + " +- " + dec.format(sdBranchs));
        dPlot.setLimits(0, histoBranchs[0].getMaximum(), 0, histoBranchs[1].getMaximum()+10);
        dPlot.draw();
        dPlot.addToStack();
        
        // Plot nbEndPoints Numbers
        ArrayUtil endPoints = new ArrayUtil(nbEndPoints);
        double meanEndPoints = endPoints.getMean();
        double sdEndPoints = endPoints.getStdDev();
        ArrayUtil[] histoEndPoints = endPoints.getHistogram(100);
        dPlot.setColor(Color.red, Color.red);
        dPlot.add("bar", histoEndPoints[0].getArray(), histoEndPoints[1].getArray());
        dPlot.addLabel(0.4, 0.2, "mean end points number = "+ dec.format(meanEndPoints) + " +- " + dec.format(sdEndPoints));
        dPlot.setLimits(0, histoEndPoints[0].getMaximum(), 0, histoEndPoints[1].getMaximum()+10);
        dPlot.draw();
        dPlot.addToStack();
        
        // Plot junctions Numbers
        ArrayUtil nbJunctions = new ArrayUtil(junctions);
        double meanJunctions = nbJunctions.getMean();
        double sdJunctions = nbJunctions.getStdDev();
        ArrayUtil[] histoJunctions = nbJunctions.getHistogram(100);
        dPlot.setColor(Color.green, Color.green);
        dPlot.add("bar", histoJunctions[0].getArray(), histoJunctions[1].getArray());
        dPlot.addLabel(0.4, 0.20, "mean junctions number = "+ dec.format(meanJunctions) + " +- " + dec.format(sdJunctions));
        dPlot.setLimits(0, histoJunctions[0].getMaximum(), 0, histoJunctions[1].getMaximum()+10);
        dPlot.draw();
        dPlot.addToStack();
        ImagePlus imgPlotDiameter = dPlot.getImagePlus();
        FileSaver plotSave = new FileSaver(imgPlotDiameter);
        plotSave.saveAsTiff(outDirResults + rootName + "_Roi"+ (r+1) + "_Plot" + ".tif");
        ResultsTable rt = dPlot.getResultsTable();
        rt.setHeading(0,"Diameters");
        rt.setHeading(1, "Branchs");
        rt.setHeading(2, "End Points");
        rt.setHeading(2, "Junctions");
        rt.updateResults();
        return(rt);
    }

    
    
    
    /**
     * 
     * @param arg
     */
    @Override
    public void run(String arg) {
        

        try {
            if (canceled) {
                IJ.showMessage(" Pluging canceled");
                return;
            }
            imageDir = IJ.getDirectory("Choose directory containing image files...");
            if (imageDir == null) {
                return;
            }
            File inDir = new File(imageDir);
            String[] imageFile = inDir.list();
            if (imageFile == null) {
                return;
            }
            // create output folder
            outDirResults = inDir + File.separator+ "Out"+ File.separator;
            File outDir = new File(outDirResults);
            if (!Files.exists(Paths.get(outDirResults))) {
                outDir.mkdir();
            }

            
            // create OME-XML metadata store of the latest schema version
            ServiceFactory factory;
            factory = new ServiceFactory();
            OMEXMLService service = factory.getInstance(OMEXMLService.class);
            IMetadata meta = service.createOMEXMLMetadata();
            ImageProcessorReader reader = new ImageProcessorReader();
            reader.setMetadataStore(meta);
            Arrays.sort(imageFile);
            int imageNum = 0;
            ArrayList<String> ch = new ArrayList();
            for (int i = 0; i < imageFile.length; i++) {
                // For all tif files
                String fileExt = FilenameUtils.getExtension(imageFile[i]);
                if ("tif".equals(fileExt)) {
                    String imageName = imageDir + imageFile[i];
                    String rootName = FilenameUtils.getBaseName(imageFile[i]);
                    // Find ROI file
                    String roi_file = imageDir+rootName+".zip";
                    if (!new File(roi_file).exists()) {
                        IJ.showMessage("No ROI file found !\nRoi file should be named as <nd or ics filename>.zip");
                        return;
                       }
                    else {
                        reader.setId(imageName);
                        int series = 0;
                        reader.setSeries(series);
                        int sizeC = reader.getSizeC();
                        imageNum++;
                        boolean showCal = false;
                        String[] channels = new String[sizeC+1];
                        if (imageNum == 1) {
                            // Check calibration
                            cal.pixelWidth = meta.getPixelsPhysicalSizeX(series).value().doubleValue();
                            cal.pixelHeight = cal.pixelWidth;
                            // problem to read calibration with nd files
                            if ((meta.getPixelsPhysicalSizeZ(series) == null) || (cal.pixelWidth == 1))
                                showCal = true;
                            else
                                cal.pixelDepth = meta.getPixelsPhysicalSizeZ(series).value().doubleValue();
                            cal.setUnit("microns");
                            System.out.println("x/y cal = " +cal.pixelWidth+", z cal = " + cal.pixelDepth);
                        }
                        
                        // find rois
                        RoiManager rm = new RoiManager(false);
                        rm.runCommand("Open", roi_file);
                        int index = 0;
                        
                        ImporterOptions options = new ImporterOptions();
                        options.setColorMode(ImporterOptions.COLOR_MODE_GRAYSCALE);
                        options.setId(imageName);
                        options.setSplitChannels(false);
                        reader.setSeries(0); 
                        options.setQuiet(true);
                        
                        /**
                         * Open channel
                         */
                        
                        // Astrocyte channel
                        System.out.println("Opening Astrocyte channel : " +rootName+".tif");
                        ImagePlus imgAstro = IJ.openImage(imageDir + imageFile[i]);
                        imgAstro.setCalibration(cal);
                        ResultsTable rtfinal = new ResultsTable();
                        // for each roi open image and crop
                        for (int r = 0; r < rm.getCount(); r++) {
                            index++;                            
                            // Find in roi name the desired top and bottom stack 
                            // roi name should be roi_number-ztop-zbottom
                            String[] regExp = rm.getName(r).split("-");
                            if (Integer.parseInt(regExp[1]) < 1)
                                regExp[1] = "1";
                            if (Integer.parseInt(regExp[2]) > imgAstro.getNSlices())
                                regExp[2] = Integer.toString(imgAstro.getNSlices());
                            int zStart = Integer.parseInt(regExp[1]);
                            int zStop = Integer.parseInt(regExp[2]);
                            
                            rm.select(imgAstro,r);
                            imgAstro.updateAndDraw();
                            Roi roiAstro = imgAstro.getRoi();
                             
                            // make substack
                            ImagePlus imgAstroCrop = new Duplicator().run(imgAstro, zStart, zStop);
                            // bug in roi manager recenter roi
                            imgAstroCrop.setRoi(roiAstro);
                            roiAstro.setLocation(0, 0);
                            imgAstroCrop.updateAndDraw();
                            roiAstro = imgAstroCrop.getRoi();
                            
                            // find background
                            double bg = find_background(imgAstroCrop);
                            
                            // get astro mask
                            ImagePlus imgAstroMask = findMask(imgAstroCrop, roiAstro, bg);
                                    
                            // get skeleton
                            ImagePlus imgAstroSkel = imgAstroMask.duplicate();
                            IJ.run(imgAstroSkel, "Skeletonize (2D/3D)", "");
                            
                            // get distance map
                            ImagePlus imgAstroMap = localThickness3D(imgAstroMask);
                            
                            // find diameter processses
                            ResultsTable rt = diameterAnalysis(imgAstroSkel, imgAstroMap, rootName, r);
                            // add results to final table
                            for (int j = 0; j < rt.getCounter(); j++) {
                                rtfinal.addValue(0, rt.getValueAsDouble(0, j));
                                rtfinal.addValue(1, rt.getValueAsDouble(1, j));
                                rtfinal.addValue(2, rt.getValueAsDouble(2, j));
                                rtfinal.addValue(3, rt.getValueAsDouble(3, j));
                            }   
                            // save map
                            IJ.run(imgAstroMap, "Enhance Contrast...", "saturated=0.3");
                            IJ.run(imgAstroMap, "Fire", "");
                            IJ.resetMinAndMax(imgAstroMap);
                            IJ.run(imgAstroMap, "Calibration Bar...", "location=[Upper Left] fill=White label=Black number=5 decimal=2 font=12 zoom=1 overlay show");
                            FileSaver astroMapFile = new FileSaver(imgAstroMap);
                            astroMapFile.saveAsTiff(outDirResults + rootName + "_Roi"+ (r+1) + "_Map.tif");
                            
                            flush_close(imgAstroCrop);
                            flush_close(imgAstroMap);
                            flush_close(imgAstroSkel);
                        }
                        rtfinal.saveAs(outDirResults + rootName + "_Results.xls");
                        rtfinal.reset();
                        flush_close(imgAstro);
                    }
                }
            }
        }
        catch (IOException | FormatException | DependencyException | ServiceException ex) {
            Logger.getLogger(Astro_FineProcesses.class.getName()).log(Level.SEVERE, null, ex);
        }
        IJ.showStatus("Process done ...");
    }
}
