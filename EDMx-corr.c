#include "EDMx-corr.h"
#include <analysis.h>
#include <formatio.h>
#include <ansi_c.h>
#include <userint.h>
#include "EDMx-corr.h"

#define MaxLines 1000
#define MAXVAR 100

#define LOGFILE "\\\\CPT\\EDManalysis\\series_res.txt"
#define PARAMLIST "\\\\CPT\\EDManalysis\\paramlist_master.txt"

#define NoPoint -99999.0
#define Pi 3.1415926
#define CorDatSuffix ".cod"
#define TRUE 1
#define FALSE 0
#define COMMENTFLAG '#'

#define REV 3
#define TIMERSTOP 2
#define LOGENTRIES 13 

char serieslog[200][150];
char prefix[20],avedir[100],corrdir[100];
double **Average_Data,**Average_Error,**Corr_Data,**Corr_Error;
int **Cut_Data,**Cut_Corr,**PlotCut,startnums[500];
int pump_flag=FALSE,fitind,loadave_flag=0,Rem,zoom_flag=FALSE;
int fit_flag=FALSE,SumAveLines,SumCodLines,PlotLines;
int bval,chanval,modeval,start_index=0,aveind=0,codind=0;
int autograph=0,autofit=0;

double *XFit_Array,*YFit_Array,xmin,xmax;
double *XFitArray=0,*YFitArray=0,**PlotError,**PlotArray;
int	   xparam,yparam,index,numchan;
double slope,intercept;

int xtype,ytype;
double x,y;
static int numlistitems=0,timercount,cursindex=0,codavenum,overlap;
static int prevcheck[2][REV]={{0,0},{0,0}}; //first index: x/y, second index place in listbox
static int itscod=0;  // plotarray uses corr scale
static int prev_xtype=-1,prev_ytype=-1,prev_xparam=-1,prev_yparam=-1,prev_Rem=0;
static double xplotstart,xplotstop,yplotstart,yplotstop;

int big,array_AveLines[1000];
int ColorArray[8] = {VAL_GREEN, VAL_YELLOW,VAL_CYAN, VAL_WHITE, VAL_BLUE, VAL_MAGENTA, VAL_GRAY, VAL_RED};

int Fill_XYArrays(int avenum); 
void Swap(double *val1, double *val2);
int Plot_points(double xstart,double xstop,double ystart,double ystop);
int FileRead(char *name, int skip, int *cols, int *rows, double **array);
void Swap(double *val1, double *val2);
int EvalKeypress(int keycode,int i);
int LoadAverage(int filenum);
int LinearFit(int lines, int maxlines, double *m, double *b, double *mdev, double *bdev, double *xsq);
double pow(double, double);
double CalcChiSqd(double *xarray,double *yarray,double *errarray,double slope,double intercept,int lines);
int SignReverse (void);
int LoadCorrelation(int filenum,int avenum);
int Cut_File_Points(void);
void MakeSpace(int totalruns,int avenum);
int compareint(const void *vp,const void *vq);

// flex param handling
char **paramlist,**readlist;
int InitParamList(int *num);
int GetFileParam(char *filename,int linenum,char *list[20]);
int MatchParam(char *name,char **list,int starti,int listnum);


// functions from linfit.c
void fit(double x[], double y[], int ndata, double sig[], double *a,
			double *b, double *siga, double *sigb, double *chi2, double *q);
void fitexy(double x[], double y[], int ndat, double sigx[], double sigy[],
			double *a, double *b, double *siga, double *sigb, double *chi2, double *q);


int main (int argc, char *argv[])
{
	int i;

	if (InitCVIRTE (0, argv, 0) == 0)	/* Needed if linking in external compiler; harmless otherwise */
		return -1;	/* out of memory */
	if ((big = LoadPanel (0, "EDMx-corr.uir", AP)) < 0)
		return -1;
	
	// flex_param
	paramlist= malloc (MAXVAR*sizeof(char *));
	for (i=0;i<MAXVAR;i++) paramlist[i] = malloc(20*sizeof(char));
	InitParamList(&numchan);
	readlist= malloc (numchan*sizeof(char *));
	for (i=0;i<numchan;i++) readlist[i] = malloc(20*sizeof(char));
	/************/
	sprintf(paramlist[numchan],"Common");
	
	Average_Data = calloc (numchan, sizeof (int *));
	for (i=0;i<numchan;++i) Average_Data[i] = calloc(10,sizeof (double));
	Average_Error = calloc (numchan, sizeof (int *));
	for (i=0;i<numchan;++i) Average_Error[i] = calloc(10,sizeof (double));
	Cut_Data = calloc (numchan, sizeof (int *));
	for (i=0;i<numchan;++i) Cut_Data[i] = calloc(10,sizeof (int));
	Corr_Data = calloc (numchan, sizeof (int *));
	for (i=0;i<numchan;++i) Corr_Data[i] = calloc(10,sizeof (double));
	Corr_Error = calloc (numchan, sizeof (int *));
	for (i=0;i<numchan;++i) Corr_Error[i] = calloc(10,sizeof (double));
	Cut_Corr = calloc (numchan, sizeof (int *));
	for (i=0;i<numchan;++i) Cut_Corr[i] = calloc(10,sizeof (int));
	XFit_Array = calloc (10,sizeof (double));
	YFit_Array = calloc (10, sizeof (double));
	PlotArray = calloc(2,sizeof(int *));
	for(i=0;i<2;++i) PlotArray[i] = calloc (10,sizeof (double));
	PlotError = calloc(2,sizeof(int *));
	for(i=0;i<2;++i) PlotError[i] = calloc (10,sizeof (double));
	PlotCut = calloc(2,sizeof(int *));
	for(i=0;i<2;++i) PlotCut[i] = calloc (10,sizeof (double));
	
	DisplayPanel (big);
	setlist(big,AP_loadlist,EVENT_COMMIT,0,0,0);
	RunUserInterface ();
	
	return 0;
}


int CVICALLBACK setlist (int panel, int control, int event,
		void *callbackData, int eventData1, int eventData2)
{
	char listname[100],line[500],label[10];
	FILE *fp;
	long size;
	int start,i,checked,filenum;
	
	if(event==EVENT_COMMIT) {
	switch (control) {
		case AP_loadlist:
		ClearListCtrl (big, AP_file_listbox);
		numlistitems=0;
		GetCtrlVal(big,AP_filelist,listname);
		if(GetFileInfo(listname,&size)) fp=fopen(listname,"r");
		else {
			ResetTextBox(big,AP_Error,"list file not found");
			return 1;
			}
		while(fgets(line,500,fp)) if(line[0] != '#') {
			sscanf(line,"%s",label);
			start=atoi(label);
			InsertListItem (big, AP_file_listbox, -1, label, start);
			CheckListItem (big, AP_file_listbox, numlistitems, 1);
			++numlistitems;
			}
		select(big,AP_B_ring,EVENT_COMMIT,0,0,0);
		select(big,AP_chan_ring,EVENT_COMMIT,0,0,0);		
		select(big,AP_corr_ring,EVENT_COMMIT,0,0,0);
		select(big,AP_HVS,EVENT_COMMIT,0,0,0);
		select(big,AP_CellNum,EVENT_COMMIT,0,0,0);
		break;
		
		case AP_savelist:
		GetCtrlVal(big,AP_filelist,listname);
		if(CompareStrings(listname,0,LOGFILE,0,0)==0) {
			MessagePopup("","Can't overwrite log file");
			return 1;
			}
		if(GetFileInfo(listname,&size)) 
			if(!ConfirmPopup("","File exists, overwrite?")) return 1;
		fp=fopen(listname,"w");
		for(i=0;i<numlistitems;++i) {
			IsListItemChecked(big,AP_file_listbox,i,&checked);
			if(checked) {
				GetLabelFromIndex(big,AP_file_listbox,i,label);
				fprintf(fp,"%s\n",label);
				}
			}
		fclose(fp);
		break;
		
		case AP_ADDTOLIST:
		GetCtrlVal(big,AP_FILENUM,&filenum);
		sprintf(label,"%d",filenum);
		InsertListItem (big, AP_file_listbox, -1, label, filenum);
		++numlistitems;
		}
	}
	return 0;
}


int CVICALLBACK select (int panel, int control, int event,
		void *callbackData, int eventData1, int eventData2)
{
	int val,i,checked,start,a[10],test,col;
	char line[500];
	float f[2];
	FILE *fp;
	
	if(event==EVENT_COMMIT) {
	GetCtrlVal(big,control,&val);

	if(val) {
		fp=fopen(LOGFILE,"r");
		for(i=0;i<numlistitems;++i) {
			IsListItemChecked(big,AP_file_listbox,i,&checked);
			if(checked) {
				GetValueFromIndex(big,AP_file_listbox,i,&start);
				do {fgets(line,500,fp);
					test=sscanf(line,"%d%d%d%d%d%d%d%f%f%d%d%d",&a[0],&a[1],&a[2],&a[3],&a[4],
								&a[5],&a[6],&f[0],&f[1],&a[7],&a[8],&a[9]);
					} while(a[0]<start || line[0]=='#');
				if(a[0] != start) {   //run# not in log file
					sprintf(line,"%d not in log",start);
					ResetTextBox(big,AP_Error,line);
					if(i<(numlistitems)) {
						GetValueFromIndex(big,AP_file_listbox,i+1,&start);
						if(a[0]==start) ++i;
						}
					}
				switch(control) {
					case AP_B_ring: col=3;
						break;
					case AP_chan_ring: col=4;
						break;
					case AP_corr_ring: col=2;
						break;
					case AP_HVS: col=7;
						if(val==0) return 1;   //0 = any switch time
						break;
					case AP_CellNum: col=8;
						if(val==0) return 1;	// 0 = any cell #
						break;
					}
				if(a[col] != val) {
					CheckListItem(big,AP_file_listbox,i,0); //uncheck
					if(control==AP_CellNum && a[col+1]==val)
						CheckListItem(big,AP_file_listbox,i,1);  // cell can be either top or bottom
					}
				}
			}
		fclose(fp);		
		}
	}	
	return 0;
}


int Plot_points(double xstart,double xstop,double ystart,double ystop)
{   
	int j,maxind,minind,show,err=0;
	double max, min;
	
	GetCtrlVal(big,AP_showcut,&show);
	DeleteGraphPlot (big, AP_Grid, -1,VAL_DELAYED_DRAW);
    SetCtrlAttribute (big, AP_Grid, ATTR_REFRESH_GRAPH, 0);
	ResetTextBox(big,AP_Error,"");
	for (j=0; j < PlotLines && err>=0; ++j) {
 		if(xstart<=PlotArray[0][j] && xstop>=PlotArray[0][j] && ystart<=PlotArray[1][j] && ystop>=PlotArray[1][j]) {
 			if (show && (PlotCut[0][j]==1 || PlotCut[1][j]==1) && PlotCut[0][j]!=2 && PlotCut[1][j]!=2) {
				err=PlotPoint(big, AP_Grid, PlotArray[0][j], PlotArray[1][j], VAL_EMPTY_SQUARE, VAL_RED);	/*plot cut points in red*/
				if(PlotError[1][j]!=0) err=PlotLine(big, AP_Grid, PlotArray[0][j],PlotArray[1][j]-PlotError[1][j],
													PlotArray[0][j], PlotArray[1][j]+PlotError[1][j], VAL_RED);
				if(PlotError[0][j]!=0) err=PlotLine(big, AP_Grid, PlotArray[0][j]-PlotError[0][j],PlotArray[1][j],
													PlotArray[0][j]+PlotError[0][j],PlotArray[1][j], VAL_RED);
				}			
   			else if(PlotCut[0][j]==0&&PlotCut[1][j]==0) {
   				err=PlotPoint(big, AP_Grid, PlotArray[0][j], PlotArray[1][j], VAL_EMPTY_SQUARE, VAL_YELLOW); /*plot uncut points in yellow*/
				if(PlotError[1][j]!=0) err=PlotLine(big, AP_Grid, PlotArray[0][j],PlotArray[1][j]-PlotError[1][j],
													PlotArray[0][j], PlotArray[1][j]+PlotError[1][j], VAL_YELLOW);
				if(PlotError[0][j]!=0) err=PlotLine(big, AP_Grid, PlotArray[0][j]-PlotError[0][j],PlotArray[1][j],
													PlotArray[0][j]+PlotError[0][j],PlotArray[1][j], VAL_YELLOW);
				}
			}
		}
	if(err==-12)
		ResetTextBox(big,AP_Error,"plot out of memory");
	if(fit_flag==TRUE) {		 /*Draw linear fit*/
		MaxMin1D(XFitArray, fitind, &max, &maxind, &min, &minind);
		PlotLine(big, AP_Grid, XFitArray[minind], slope*XFitArray[minind]+intercept, XFitArray[maxind],
						slope*XFitArray[maxind]+intercept, VAL_RED);
		}
	SetCtrlAttribute (big,AP_Grid, ATTR_REFRESH_GRAPH, 1);
	
	return 0;
}       		    


int CVICALLBACK PlotSet (int panel, int control, int event, void *callbackData, int eventData1, int eventData2)
{						/*the "command and control" function*/
	int imin,newlen,i,j,select,scan;
	static int mindind,maxdind,sequence=0,direction=0,show;
	double ymin,ymax,*Massldiff,xzoommin,xzoommax,yzoommin,yzoommax;
	static double mind,maxd;
	char graphlabel[100],filename[50],filenum[20];


	GetCtrlVal(big,AP_XParameter,&xparam);
	GetCtrlVal(big,AP_YParameter,&yparam);
	GetCtrlVal(big,AP_Removal,&select);
	GetCtrlVal(big,AP_XCorrAve,&xtype);
	GetCtrlVal(big,AP_YCorrAve,&ytype);
	
	switch (event) {
		case EVENT_COMMIT:
			if(control==AP_Graph_button) {	/*graph button pressed*/				
				fit_flag=FALSE;
				GetCtrlVal(big,AP_remove,&Rem);
				GetCtrlVal(big,AP_showcut,&show);
       			if(Rem!=prev_Rem) Cut_File_Points();
     			prev_Rem=Rem;
       			Fill_XYArrays(codavenum);
				SignReverse();
				MaxMin1D(PlotArray[0],PlotLines,&xplotstop,&maxdind,&xplotstart,&mindind);
				MaxMin1D(PlotArray[1],PlotLines,&yplotstop,&maxdind,&yplotstart,&mindind);
       			Plot_points(xplotstart,xplotstop,yplotstart,yplotstop);
       			}
			else if(control==AP_Reset) {		 /*reset button pressed - view plot over entire area*/
				fit_flag=FALSE;
				MaxMin1D(PlotArray[0],PlotLines,&xplotstop,&maxdind,&xplotstart,&mindind);
				MaxMin1D(PlotArray[1],PlotLines,&yplotstop,&maxdind,&yplotstart,&mindind);
       			Plot_points(xplotstart,xplotstop,yplotstart,yplotstop);
       			}
       		else if(control==AP_Zoom) {			   /*zoom button pressed - scope in*/
       			fit_flag=FALSE;
       			GetGraphCursor (big, AP_Grid, 1, &xplotstart, &yplotstart);
				GetGraphCursor (big, AP_Grid, 2, &xplotstop, &yplotstop);
				if(xplotstop<xplotstart) Swap(&xplotstart,&xplotstop);
				if(yplotstop<yplotstart) Swap(&yplotstart,&yplotstop);
       			Plot_points(xplotstart,xplotstop,yplotstart,yplotstop);
       			}
       		break;									 

       	case EVENT_LEFT_CLICK:
       		/* in single cursor mode: activate timer to get cursor position */
       		/* after cursor snaps to point									*/
       		if(control==AP_Grid && select==1) {
       			SetCtrlAttribute(big,AP_CursorTimer,ATTR_ENABLED,1);
       			sequence=0;
       			timercount=0;
       			}
			break;
		
		case EVENT_KEYPRESS:					  /*clear,add,or cut point in cursor mode*/
       		if(control==AP_Grid && select==1) {
				timercount=TIMERSTOP;
				if(!sequence) {
					GetGraphCursor (big, AP_Grid, 1, &x, &y);
					Massldiff=calloc(PlotLines,sizeof(double));
					for(i=0;i<PlotLines;i++) Massldiff[i]=fabs(PlotArray[0][i]-x);
					MaxMin1D(Massldiff,PlotLines,&maxd,&maxdind,&mind,&cursindex);
					fit_flag=FALSE;
					sequence=1;
					free(Massldiff);
					}
				EvalKeypress(eventData1,cursindex);			
				Plot_points(xplotstart,xplotstop,yplotstart,yplotstop);
				if(eventData1==88 || eventData1==65 || eventData1==2304) { //"Shift"=forward
					direction=1;
					}
				else if(eventData1==262232 || eventData1==262209 || eventData1==2048) {//"Ctrl"=backward
					direction=-1;
					}
				else direction=0;
				cursindex += direction;
				if(direction)
					if(!show) while(PlotCut[0][cursindex] || PlotCut[1][cursindex]) {
						cursindex += direction;
						}
					else while(PlotCut[0][cursindex]==2 || PlotCut[1][cursindex]==2) {
						cursindex += direction;
						}
				x = PlotArray[0][cursindex]; y=PlotArray[1][cursindex]; scan=Corr_Data[0][cursindex];
				SetGraphCursor (big, AP_Grid,1,x,y);
				sprintf(graphlabel,"x = %g y = %g  scan# %d",x,y,scan);
				SetCtrlAttribute(big, AP_Grid,ATTR_LABEL_TEXT,graphlabel);
				}
			else if(control==AP_Grid && select==2) {	   /*clear,add,or cut points in area mode*/
				fit_flag=FALSE;
				GetGraphCursor (big, AP_Grid, 1, &xmin, &ymin);
				GetGraphCursor (big, AP_Grid, 2, &xmax, &ymax);
				if(xmin>xmax) Swap(&xmin,&xmax);
				if(ymin>ymax) Swap(&ymin,&ymax);
				for(j=0;j<PlotLines;++j) {
					x=PlotArray[0][j];y=PlotArray[1][j];
					if((eventData1==262232 || eventData1==262209 || eventData1==262213) && PlotCut[0][j]!=2 && PlotCut[1][j]!=2) {
						if(!(x>xmin && x<xmax) || !(y>ymin && y<ymax)) 
							EvalKeypress(eventData1,j);
						}
					else if((x>xmin && x<xmax) && (y>ymin && y<ymax)) 
						EvalKeypress(eventData1,j);
					}
				Plot_points(xplotstart,xplotstop,yplotstart,yplotstop);
				GetGraphCursor (big, AP_Grid, 1, &xmin, &ymin);
				GetGraphCursor (big, AP_Grid, 2, &xmax, &ymax);
				if(autofit) Draw_Fit(big,AP_LinearFit,EVENT_COMMIT,0,0,0);
				}
			break;
		}
	return 0;
}


void Swap(double *val1, double *val2)		 /*swap max, min values if max<min*/
{	
	double valsave;
	
	valsave = *val1; *val1=*val2; *val2=valsave;
}


int CVICALLBACK Quit (int panel, int control, int event,
		void *callbackData, int eventData1, int eventData2)
{
	if(event==EVENT_COMMIT) {    /*free memory, close program*/
		QuitUserInterface (0);
		free(Average_Data);
		free(Cut_Data);
		free(Cut_Corr);
		free(Corr_Data);
		free(Corr_Error);
		free(Average_Error);
		free(PlotArray);
		free(PlotError);
		free(PlotCut);
		free(XFit_Array);
		free(YFit_Array);
		}
	return 0;
}


int FileRead(char *name, int skip, int *cols, int *rows,double **array)
{  
	char line[3000], *tmp;
	int i,j;
	FILE *fp;
	fpos_t ptr;
	
	ProcessSystemEvents;
	if (!(fp=fopen(name,"r"))) return 1;
	*rows=0;
	while (fgets(line,3000,fp)) (*rows)++;
	rewind(fp);    
	for (i=0;i<skip;i++) fgets(line,3000,fp);
	*rows=*rows-skip;
	fgetpos(fp,&ptr);
	fgets(line,3000,fp);
	strtok(line," \t,");
	*cols=1;
	while (strtok(NULL," \t,")) (*cols)++;
	*array= realloc(*array,(*cols)*(*rows)*sizeof(double));
	fsetpos(fp,&ptr);
	i=0;
	while (fgets(line,3000,fp)) {
		j=0;
		(*array)[j*(*rows)+i]=atof(strtok(line," \t,"));
		while (tmp=strtok(NULL," \t,")) {
			j++;
			(*array)[j*(*rows)+i]=atof(tmp);
			}
		++i;
		}
	fclose(fp);
	
	return 0;
}


int EvalKeypress(int keycode, int q)
{
	switch(keycode) {
		case 120: case 88: case 262232:		// x or X or (Ctrl)x
			PlotCut[0][q]=1;
			PlotCut[1][q]=1;
			break;
		case 97: case 65: case 262209:  	// a or A or (Ctrl)a
			PlotCut[0][q]=0;
			PlotCut[1][q]=0;
			break;
		case 101: case 69: case 262213:		// e or E or (Ctrl)e
			PlotCut[0][q]=2;
			PlotCut[1][q]=2;
			break;
			default: return 1;	//invalid key
		}
	return 0;
}


int CVICALLBACK RemovePoints (int panel, int control, int event,
		void *callbackData, int eventData1, int eventData2)
{
	int mode,i,j;								  /*Set cursor attributes according to remove mode*/
	static double x1,x2,y1,y2;
	static int prev_xparam=0,prev_yparam=0,prev_mode=0;

	if(event==EVENT_COMMIT) {
		
		switch (control) {
		case AP_Removal:
			GetCtrlVal(big,AP_Removal,&mode);
			if(prev_mode!=1) {
				GetGraphCursor (big, AP_Grid, 1, &x1, &y1);
				GetGraphCursor (big, AP_Grid, 2, &x2, &y2);
				}
			if(mode!=1) {
				SetCtrlAttribute (big, AP_Grid, ATTR_LABEL_TEXT, "");
				SetCtrlAttribute (big, AP_Grid, ATTR_NUM_CURSORS, 2);
				for(i=1;i<=2;++i){
					SetCursorAttribute (big, AP_Grid, i, ATTR_CROSS_HAIR_STYLE,
										VAL_LONG_CROSS);
					SetCursorAttribute (big, AP_Grid, i, ATTR_CURSOR_MODE,
										VAL_FREE_FORM);
					SetCursorAttribute (big, AP_Grid, i, ATTR_CURSOR_POINT_STYLE,
									VAL_SMALL_EMPTY_SQUARE);
					SetCursorAttribute (big, AP_Grid, i, ATTR_CURSOR_COLOR,
										ColorArray[3+i]);
					}
				if(prev_xparam==xparam && prev_yparam==yparam) {
					SetGraphCursor (big, AP_Grid, 1, x1, y1);
					SetGraphCursor (big, AP_Grid, 2, x2, y2);
					}
				SetCtrlAttribute (big, AP_Zoom, ATTR_DIMMED, 0);
				}
			else {
				SetCtrlAttribute (big, AP_Grid, ATTR_NUM_CURSORS, 1);
				SetCursorAttribute (big, AP_Grid, 1, ATTR_CURSOR_COLOR, VAL_RED);
				SetCursorAttribute (big, AP_Grid, 1, ATTR_CROSS_HAIR_STYLE,
									VAL_SHORT_CROSS);
				SetCursorAttribute (big, AP_Grid, 1, ATTR_CURSOR_MODE,
									VAL_SNAP_TO_POINT);
				SetCtrlAttribute (big, AP_Zoom, ATTR_DIMMED, 1);
				}
			prev_mode=mode; prev_xparam=xparam; prev_yparam=yparam;
			break;
			}
	}
	return 0;
}


int CVICALLBACK Draw_Fit (int panel, int control, int event,
		void *callbackData, int eventData1, int eventData2)
{
	int i,avenum,CurrentCtrl;										   
	static double m, b, mdev, bdev, chi2,q;
	double errmult=1;
	
	if(event==EVENT_COMMIT) {
		fit_flag=TRUE;
		XFitArray = realloc (XFitArray,2*PlotLines*sizeof(double));
		YFitArray = realloc (YFitArray,2*PlotLines*sizeof(double));
		fitind=0;
		for(i=0;i<PlotLines;++i) {
			if(PlotCut[0][i]==0 && PlotCut[1][i]==0) {
				XFitArray[fitind]=PlotArray[0][i];
				XFitArray[PlotLines+fitind]=PlotError[0][i];
				YFitArray[fitind]=PlotArray[1][i];
				YFitArray[PlotLines+fitind]=PlotError[1][i];
				++fitind;
				}
			}
		if (fitind==0) return -1;
		LinearFit (fitind, PlotLines, &m, &b, &mdev, &bdev, &chi2);
		slope=m; intercept=b;
		Plot_points(xplotstart,xplotstop,yplotstart,yplotstop);
		if(overlap) {
			GetCtrlVal(big,AP_NumToAve,&avenum);
			switch(avenum) {	 //correction for non-independent errors
				case 3: errmult = 4/sqrt(6); break;
				case 4: errmult = 8/sqrt(20); break;
				case 5: errmult = 16/sqrt(70); break;
				case 6: errmult = 32/sqrt(252); break;
				}
			mdev *= errmult;
			bdev *= errmult;
			}
		
		SetCtrlVal(big,AP_Slope,m);
		SetCtrlVal(big,AP_Intercept,b);
		SetCtrlVal(big,AP_Slope_Error,mdev);
		SetCtrlVal(big,AP_Int_Error,bdev);
		SetCtrlVal(big,AP_SlopeOverError,m/mdev);
		if(fabs(m/mdev) > 2.0) SetCtrlAttribute(big,AP_SlopeOverError,ATTR_FRAME_COLOR,VAL_RED);
		else SetCtrlAttribute(big,AP_SlopeOverError,ATTR_FRAME_COLOR,VAL_LT_GRAY);
		SetCtrlVal(big,AP_IntOverError,b/bdev);
		if(fabs(b/bdev) > 2.0) SetCtrlAttribute(big,AP_IntOverError,ATTR_FRAME_COLOR,VAL_RED);
		else SetCtrlAttribute(big,AP_IntOverError,ATTR_FRAME_COLOR,VAL_LT_GRAY);
		SetCtrlVal(big,AP_chi2,chi2);
		CurrentCtrl = GetActiveCtrl(big);	// enables quick scrolling through parameters when in autograph mode
											// need to return to the main panel and last control after stdin is opened
		//other fit methods to check xy error weighting	
		// linear fit using method from Numerical Recipes, first with x and y error bars, then with only y errors
		fitexy(XFitArray,YFitArray,fitind,&XFitArray[PlotLines],&YFitArray[PlotLines],&b,&m,&bdev,&mdev,&chi2,&q);
		chi2=chi2/(fitind-2);
		printf("%13s%13s%13s%13s%13s%13s\n","slope","err","intercept","err","chi2","q");
		printf("%13g%13g%13g%13g%13g%13g\n",m,mdev*sqrt(chi2)*errmult,b,bdev*sqrt(chi2)*errmult,chi2,q);
		fit(XFitArray,YFitArray,fitind,&YFitArray[PlotLines],&b,&m,&bdev,&mdev,&chi2,&q);
		chi2=chi2/(fitind-2);
		printf("%13g%13g%13g%13g%13g%13g\n*****************\n",m,mdev*sqrt(chi2)*errmult,b,bdev*sqrt(chi2)*errmult,chi2,q);		
		
		SetActivePanel(big);
		SetActiveCtrl(big,CurrentCtrl);
		}
	return 0;
}


int LinearFit(int lines, int maxlines, double *m, double *b, double *mdev, double *bdev,double *xsq)
{	/*calculate linear fit parameters*/
	double delta,s1,s2,s3,s4,s5,dum=0,*weighted_dev;
	int i,j,weightchisq;

	weighted_dev=calloc(maxlines,sizeof(double));
	memset(weighted_dev,0,maxlines*sizeof(double));

	if(YFitArray[maxlines]==0.0 && XFitArray[maxlines]!=0.0) {   // x error bars, no y errors: invert x and y
		s1=0; s2=0; s3=0; s4=0; s5=0;
		for(i=0;i<lines;++i) {
			s1 += pow(YFitArray[i],2)/pow(XFitArray[maxlines+i],2);
			s2 += (XFitArray[i])/pow(XFitArray[maxlines+i],2);
			s3 += (YFitArray[i])/pow(XFitArray[maxlines+i],2);
			s4 += (YFitArray[i])*(XFitArray[i])/pow(XFitArray[maxlines+i],2);
			s5 += 1/pow(XFitArray[maxlines+i],2);
			}
		delta=s5*s1-pow(s3,2);

		*b=(s1*s2-s3*s4)/delta;
		*m=(s5*s4-s3*s2)/delta;
		*bdev=sqrt(s1/delta);
		*mdev=sqrt(s5/delta);
		*xsq = CalcChiSqd(YFitArray,XFitArray,&XFitArray[maxlines],*m,*b,lines);
		if(*xsq>1) {
			*bdev=(*bdev)*sqrt(*xsq);
			*mdev=(*mdev)*sqrt(*xsq);
			}
		*mdev=*mdev/pow(*m,2);
		*bdev=sqrt(pow(*bdev,2)+pow(*b,2)/pow(*m,2)*pow(*mdev,2))/(*m);
		*b=(-1)*(*b)/(*m);
		*m=1/(*m);
		}
	else {
		if(YFitArray[maxlines]==0.0) {  // no error bars: equal weighting
			for(i=0;i<lines;++i) YFitArray[maxlines+i]=1.0;
			weightchisq=1;
			}
		s1=0; s2=0; s3=0; s4=0; s5=0;
		for(i=0;i<lines;++i) {
			weighted_dev[i]=YFitArray[maxlines+i];
			s1 += pow(XFitArray[i],2)/pow(YFitArray[maxlines+i],2);
			s2 += (YFitArray[i])/pow(YFitArray[maxlines+i],2);
			s3 += (XFitArray[i])/pow(YFitArray[maxlines+i],2);
			s4 += (XFitArray[i])*(YFitArray[i])/pow(YFitArray[maxlines+i],2);
			s5 += 1/pow(YFitArray[maxlines+i],2);
			}
		delta=s5*s1-pow(s3,2);
		*b=(s1*s2-s3*s4)/delta;
		*m=(s5*s4-s3*s2)/delta;
	
		if(XFitArray[maxlines]!=0.0) for(j=0;j<20;++j) {   // x and y error bars: propagate x errors into weighting with fitted slope
			s1=0; s2=0; s3=0; s4=0; s5=0;
			for(i=0;i<lines;++i) {
				weighted_dev[i]=sqrt(pow(YFitArray[maxlines+i],2) + pow(XFitArray[maxlines+i]*(*m),2));
				s1 += pow(XFitArray[i],2)/pow(weighted_dev[i],2);
				s2 += (YFitArray[i])/pow(weighted_dev[i],2);
				s3 += (XFitArray[i])/pow(weighted_dev[i],2);
				s4 += (XFitArray[i])*(YFitArray[i])/pow(weighted_dev[i],2);
				s5 += 1/pow(weighted_dev[i],2);
				}
			delta=s5*s1-pow(s3,2);
			*b=(s1*s2-s3*s4)/delta;
			*m=(s5*s4-s3*s2)/delta;
			}
		*bdev=sqrt(s1/delta);
		*mdev=sqrt(s5/delta);
		*xsq=CalcChiSqd(XFitArray,YFitArray,weighted_dev,*m,*b,lines);
		if(*xsq>1 || weightchisq) {
			*bdev=*bdev*sqrt(*xsq);
			*mdev=*mdev*sqrt(*xsq);
			}
		}
	return 0;
}


double CalcChiSqd(double *xarray,double *yarray,double *errarray,double slope,double intercept,int lines)
{
	double fit,xsq, dum=0;
	int i;

	for(i=0;i<lines;++i){
		fit=slope*xarray[i]+intercept;
		dum=dum+pow((yarray[i]-fit)/errarray[i],2);
		}
	xsq=dum/(lines-2);

	return xsq;
}


int CVICALLBACK CheckAllRemoveAll (int panel, int control, int event,
		void *callbackData, int eventData1, int eventData2)
{								/* choose all or no files */
	int tot_count,i;

	if(event==EVENT_COMMIT) {
		GetNumListItems (big, AP_file_listbox, &tot_count);   
		if(control==AP_checkall) {
			for(i=0;i<tot_count;++i) CheckListItem (big, AP_file_listbox, i, 1);
			}
		else if(control==AP_removeall) {
			for(i=0;i<tot_count;++i) CheckListItem (big, AP_file_listbox, i, 0);
			}
		}
	return 0;
}


int SignReverse(void)
{						/*change sign of data based on magnetic field, channel signs*/
	int i,j,k,n,start=0,stop=0,a[3],check[2][REV],logpos=0,topcell,botcell,test;
	double f[5];
	char line[1000],box[5],other[5][20],junk[50];

	for(i=0;i<2;++i) {		// x/y
		for(j=0;j<REV;++j) {  // item in listbox
			if(i==0) IsListItemChecked(big,AP_XSignRev,j,&check[i][j]);
			else IsListItemChecked(big,AP_YSignRev,j,&check[i][j]);
			if(prevcheck[i][j]!=check[i][j]) {
				logpos=0;
				test=sscanf(serieslog[logpos],"%d%d%d%d%d%lg%lg%lg%lg%lg%d%d%s%[ \t]%[^;\n];%[^;\n];%[^;\n];%[^;\n];%[^;\n]",&start,&stop,&a[0],&a[1],&a[2],
							&f[0],&f[1],&f[2],&f[3],&f[4],&topcell,&botcell,box,junk,other[0],other[1],other[2],other[3],other[4]);
				for(n=0; n+LOGENTRIES+1 < test; ++n) {
					if(!CompareStrings(other[n],0,"HVT",0,1)) f[3] *= 0.5;
					else if(!CompareStrings(other[n],0,"HVB",0,1)) f[3] *= 0.5;
					}
				for(k=0;k<PlotLines;++k) {
					if((itscod && Corr_Data[0][k]>stop) || (!itscod && Average_Data[0][k]>stop)) {
						++logpos;
						test=sscanf(serieslog[logpos],"%d%d%d%d%d%lg%lg%lg%lg%lg%d%d%s%[ \t]%[^;\n];%[^;\n];%[^;\n];%[^;\n];%[^;\n]",&start,&stop,&a[0],&a[1],&a[2],
									&f[0],&f[1],&f[2],&f[3],&f[4],&topcell,&botcell,box,junk,other[0],other[1],other[2],other[3],other[4]);
						for(n=0; n+LOGENTRIES+1 < test; ++n) {
							if(!CompareStrings(other[n],0,"HVT",0,1)) f[3] *= 0.5;
							else if(!CompareStrings(other[n],0,"HVB",0,1)) f[3] *= 0.5;
							}
						}
					if(j==2) {
						if(check[i][2]) {
							PlotArray[i][k] *= (10.0/f[3]);
							PlotError[i][k] *= (10.0/f[3]);
							}
						else {
							PlotArray[i][k] *= (f[3]/10.0);
							PlotError[i][k] *= (f[3]/10.0);
							}
						}
					else PlotArray[i][k] *= a[1+j];
					}
				prevcheck[i][j]=check[i][j];
				}
			}
		}
	return 0;
}

	
int CVICALLBACK checkrange (int panel, int control, int event,
		void *callbackData, int eventData1, int eventData2)
{
	int i,min,max,value,check=0;
	
	if(event==EVENT_COMMIT) {
	GetCtrlVal(big,AP_min,&min);
	GetCtrlVal(big,AP_max,&max);
	
	if(control==AP_check_range) check=1;
	else if(control==AP_remove_range) check=0;
	for(i=0;i<numlistitems;++i) {
		GetValueFromIndex(big,AP_file_listbox,i,&value);
		if(min<=value && value<=max) CheckListItem(big,AP_file_listbox,i,check);
		}
	
	if(check) {
		select(big,AP_B_ring,EVENT_COMMIT,0,0,0);
		select(big,AP_chan_ring,EVENT_COMMIT,0,0,0);		
		select(big,AP_corr_ring,EVENT_COMMIT,0,0,0);
		select(big,AP_HVS,EVENT_COMMIT,0,0,0);
		select(big,AP_CellNum,EVENT_COMMIT,0,0,0);
		}
	}
	return 0;
}


int CVICALLBACK Load_Stuff (int panel, int control, int event,		 /*load data from chosen files*/
		void *callbackData, int eventData1, int eventData2)
{
	int i,checked,dum,w,test;
	int start,stop,totalruns=0,count,cols,lines;
	char filename[100],line[150],*checkend;
	double *CumData;
	FILE *log;

	if(event==EVENT_COMMIT) {
		prev_xtype=-1,prev_ytype=-1,prev_xparam=-1,prev_yparam=-1;
		prev_Rem=0;
		GetCtrlVal(big,AP_NumToAve,&codavenum);
		GetCtrlVal(big,AP_CORRMODE,&overlap);
		GetCtrlVal(big,AP_PREFIX,prefix);
		GetCtrlVal(big,AP_AVEDIR,avedir);
		GetCtrlVal(big,AP_CORRDIR,corrdir);
		index=0;
		GetNumListItems (big, AP_file_listbox, &count);
		dum=0;
		log=fopen(LOGFILE,"r");
		for(w=0;w<count;++w) {
			IsListItemChecked (big, AP_file_listbox, w, &checked);
			if(checked==1) {
				GetValueFromIndex(big,AP_file_listbox,w,&startnums[dum]);
				do {checkend=fgets(line,150,log);
					test=sscanf(line,"%d%d",&start,&stop);
					} while((start<startnums[dum] || line[0]=='#') && checkend!=NULL);
				if(start != startnums[dum]) {   //run# not in log file
					sprintf(line,"%d not in log",startnums[dum]);
					ResetTextBox(big,AP_Error,line);
					rewind(log);
					sprintf(line,"%s%s%d-.cum",avedir,prefix,startnums[dum]);
					CumData = calloc(1,sizeof(double));
					FileRead(line,2,&cols,&lines,&CumData); //need to find number of scans
					start=CumData[0];
					stop=CumData[lines-1];
					sprintf(line,"%d %d *",start,stop);  //place dummy line in the logarray and continue
					free(CumData);
					}
				strcpy(serieslog[dum],line);   //save logline
				totalruns += (stop-start)+1;
				++dum;
				}
			}
		index=dum;
		MakeSpace(totalruns,codavenum);
		aveind=0;
		codind=0;
		for(w=0;w<index;++w) {
			SetCtrlVal(big,AP_file_being_read,startnums[w]);
			LoadAverage(startnums[w]);
			LoadCorrelation(startnums[w],codavenum);
			}
		}
	return 0;
}


int InitParamList(int *num)
// fill master parameter list and initialize single shot parameter ring
{
	int i=0,j=0,singlisti=0,type,len;
	char line[100],*name;
	FILE *fp;

	fp=fopen(PARAMLIST,"r");
	ClearListCtrl (big, AP_XParameter);
	ClearListCtrl (big, AP_YParameter);
	while (fgets(line,100,fp)) {
		if (line[0] != COMMENTFLAG && strlen(line)>18) {
			name=strtok(line,":");
			InsertListItem(big,AP_XParameter,-1,name,i);
			InsertListItem(big,AP_YParameter,-1,name,i);
			strcpy(paramlist[i],name);
			++i;
			}
		}
	*num = i;
	fclose(fp);
	
	return 0;
}


int GetFileParam(char *filename,int linenum,char *list[20])
// get first line of data file which should contain parameter names separated by at least 2 spaces
//	fills an array with the names
{
	int i=0,j,listi=0,len;
	char line[3000],tmp[18];
	FILE *fp;
	
	if (!(fp=fopen(filename,"r"))) return -1;
	for(j=0;j<linenum;++j) fgets(line,3000,fp);
	fclose(fp);
	
	len=strlen(line)-1; // length not including '\n'
	while(line[i]==' ' || line[i]==COMMENTFLAG && i<len) ++i;
	while(i< (len-3)) {  // assume names have at least 4 characters
		for (j=0; (line[j+i+1] != ' ' || line[j+i+2] != ' ') && line[j+i+2]!='\n'; j++); 
		if(j+i+2 == len) ++j; // make sure last parameter has enough space
		list[listi]=calloc(18,sizeof(char));
		strncpy(list[listi],&line[i],j+1);
		++listi;
		i += j+1;
		while((line[i]==' ' || line[i]==COMMENTFLAG) && i<len) ++i;
		}

	return listi;
}


int MatchParam(char *name,char **list,int starti,int listnum)
// find index of 'name' from list, index is returned if succesfull otherwise returns -1
{
	int i;
	
	for(i=starti;i<listnum && strcmp(name,list[i]);++i);
	if(i==listnum) for(i=0;i<starti && strcmp(name,list[i]);++i); // search from beginning

	if(strcmp(name,list[i])) return -1;  // parameter not found
	else return i; // index of name
}


void MakeSpace(int totalruns,int avenum)
{
	int i,CodLines;
	
	if(overlap) CodLines=totalruns;  // same as total for overlap
	else CodLines=totalruns/avenum;
	for(i=0;i<numchan;++i) {
		Average_Data[i] = realloc(Average_Data[i],totalruns*sizeof(double));
		Average_Error[i]=realloc(Average_Error[i],totalruns*sizeof(double));
		Cut_Data[i]=realloc(Cut_Data[i],totalruns*sizeof(int));
		memset(Cut_Data[i],0,totalruns*sizeof(int));
		Corr_Data[i]=realloc(Corr_Data[i],CodLines*sizeof(double));
		Cut_Corr[i]=realloc(Cut_Corr[i],CodLines*sizeof(int));
		memset(Cut_Corr[i],0,CodLines*sizeof(int));
		Corr_Error[i]=realloc(Corr_Error[i],CodLines*sizeof(double));
		}
	for(i=0;i<2;++i) PlotArray[i]=realloc(PlotArray[i],totalruns*sizeof(double));
	for(i=0;i<2;++i) PlotError[i]=realloc(PlotError[i],totalruns*sizeof(double));
	for(i=0;i<2;++i) PlotCut[i]=realloc(PlotCut[i],totalruns*sizeof(int));
}
		

int LoadAverage(int filenum)					 /*load average data from files chosen into array*/
{   
	int i,n,skip=0,col,AveLines,exists,dum,numparam;
	int lastp=0,masteri;
	FILE *fp;
	char line[3000],filename[100];
	double *ReadAve=0;
	
	ReadAve=malloc(sizeof(double));
	sprintf(filename,"%s%s%d%s",avedir,prefix,filenum,"-.ave");
	exists=GetFileInfo(filename,&dum);
	if(exists) {
		fp=fopen(filename,"r");
		fgets(line,3000,fp);
		while(line[0]=='#') {
			fgets(line,3000,fp);
			++skip;
			}
		fclose(fp);
		numparam=GetFileParam(filename,skip,readlist);
		FileRead(filename,skip,&col,&AveLines,&ReadAve);
		for(i=0;i<numparam;++i) {
			masteri=MatchParam(readlist[i],paramlist,lastp,numchan);
			if(masteri>=0) {
				for(n=0;n<AveLines;++n) {
					Average_Data[masteri][aveind+n]=ReadAve[2*i*AveLines+n];
					Average_Error[masteri][aveind+n]=ReadAve[2*i*AveLines+AveLines+n];
				}
				lastp=masteri;
			}
			else {
				sprintf(line,"%d ave parameter \"%s\" not in masterlist",filenum,paramlist[i]);
				ResetTextBox(big,AP_Error,line);
			}
			lastp=masteri;
		}
		aveind += AveLines;
	}
	else {
		sprintf(line,".ave file %d not found",filenum);
		ResetTextBox(big,AP_Error,line);
	}
	free(ReadAve);
		
	return 0;
}
		

int LoadCorrelation(int filenum,int avenum)				   /*load correlation data from chosen files*/
{ 
	int skip=0,col,i,j,n,CodLines,exists,dum;
	int numparam,masteri,lastp=0;
	FILE *fp;			    
	char filename[100],line[1000];
	double *ReadAve=0,sum;
			    
	ReadAve=malloc(sizeof(double));
	sprintf(filename,"%s%s%d-%d%s",corrdir,prefix,filenum,avenum,".cod");
	skip=1;
	exists=GetFileInfo(filename,&dum);
	if(exists) {
		numparam=GetFileParam(filename,skip,readlist);
		FileRead(filename,skip,&col,&CodLines,&ReadAve);
		for(i=0;i<numparam;++i) {
			masteri=MatchParam(readlist[i],paramlist,lastp,numchan);
			if(masteri>=0) {
				for(n=0;n<CodLines;++n) {
					Corr_Data[masteri][codind+n]=ReadAve[2*i*CodLines+n];
					if(Corr_Data[masteri][codind+n]==NoPoint) Cut_Corr[masteri][codind+n]=2;  // don't display
					Corr_Error[masteri][codind+n]=ReadAve[2*i*CodLines+CodLines+n];
				}
				lastp=masteri;
			}
			else {
				sprintf(line,"%d cod parameter \"%s\" not in masterlist",filenum,paramlist[i]);
				ResetTextBox(big,AP_Error,line);
			}
		}
		codind += CodLines;
	}
	else {
		sprintf(line,".cod file %d not found",filenum);
		ResetTextBox(big,AP_Error,line);
		}
	free(ReadAve);

	return 0;
}


int Cut_File_Points()			 /*include or remove cut data from file*/
{	
	int i,ind,j,k,l,startnum,stop,place=0,com=0,numcut,cut[1000];
	int exists,dum;
	char remname[55],line[MaxLines*10+20],c,lab[20],*start,*tmp;
	double *ReadAve;
	FILE *fp;
	
	if(Rem==1) {
		prev_xtype=-1,prev_ytype=-1,prev_xparam=-1,prev_yparam=-1;
		for(ind=0;ind<index;++ind) {
			sscanf(serieslog[ind],"%d %d",&startnum,&stop);
			sprintf(remname,"%s%s%d%s",avedir,prefix,startnum,"-.rem");
			exists=GetFileInfo(remname,&dum);
			if(exists) {
				fp=fopen(remname,"r");
				while (fgets(line,MaxLines*10+20,fp)) {
					if(line[0] != '#') {
						for(j=0;(c=line[j])!=':';++j);
						lab[0]=NULL; strncat(lab,line,j);
						for(j=0;CompareStrings(lab,0,paramlist[j],0,0) && j<=numchan;++j);
						if(j>numchan) ResetTextBox(big,AP_Error,"invalid parameter in .rem file");
						if(j==numchan) {  // Common
							com=1;
							j=1;
							}
					//	--j;
						start=strrchr(line,':');
						numcut=0;
						if(start!=NULL) {
							tmp=strtok(start+1," \t,");
							if(!(tmp[0]=='[')) {
								cut[numcut]=atoi(tmp);
								++numcut;
								}
							while(tmp=strtok(NULL," \t,")) {
				    		   	if(!(tmp[0]=='[')) {
				    		  		cut[numcut]=atoi(tmp);
									++numcut;
									}
				 				}
							do {
								qsort(cut,numcut,sizeof(int),compareint);
								k=0;
								while(cut[k]<startnum && k<numcut) ++k;
								for(l=0;l<(stop-startnum+1) && (k<numcut);++l) {
									if(Average_Data[0][place+l]==cut[k]) {
										Cut_Data[j][place+l]=1;
										++k;
										}
									}
								++j;
								} while(com && j<numchan);							
							}
						}
			 		}
			 	place += (stop-startnum+1);     
			 	fclose(fp);
				}
			else {
				sprintf(line,".rem file %d not found",startnum);
				ResetTextBox(big,AP_Error,line);
				}
			}
		}
	return 0;
}


int compareint(const void *vp,const void *vq)
{
	int *p=vp, *q=vq;
	return (abs(*p)-abs(*q));
}

		
int Fill_XYArrays(int avenum)
{
	int i,j,k=0,n,type[2],param[2];
	double sumave,sumerr;
	
	GetCtrlVal(big,AP_XParameter,&param[0]);
	GetCtrlVal(big,AP_YParameter,&param[1]);
	GetCtrlVal(big,AP_XCorrAve,&type[0]);
	GetCtrlVal(big,AP_YCorrAve,&type[1]);
	
	if(type[0]!=prev_xtype || type[1]!=prev_ytype || param[0]!=prev_xparam || param[1]!=prev_yparam) {
		for(i=0;i<2;++i) for(j=0;j<REV;++j) prevcheck[i][j]=0;
		if(type[0] || type[1]) {
			PlotLines=codind;
			itscod=1;
			}
		else {
			PlotLines=aveind;
			itscod=0;
			}
		for(j=0;j<2;++j) {
			for(i=0;i<PlotLines;++i) {
				if(itscod) {
					if(type[j]) {
						PlotArray[j][i]=Corr_Data[param[j]][i];
						PlotError[j][i]=Corr_Error[param[j]][i];
						PlotCut[j][i]=Cut_Corr[param[j]][i];
						}
					else {
						while(Average_Data[0][k] < Corr_Data[0][i]) ++k;
					/*	sumave=0;  // use for non overlap cod
						sumerr=0;
						for(n=0;n<avenum;++n) {
							sumave += Average_Data[param[j]][k+n];
							sumerr += Average_Error[param[j]][k+n];
							}
						PlotArray[j][i]=sumave/avenum;
						PlotError[j][i]=sumerr/avenum;	   */
						PlotArray[j][i]=Average_Data[param[j]][k];			// overlap cod
						PlotError[j][i]=Average_Error[param[j]][k];
						PlotCut[j][i]=Cut_Data[param[j]][i];
						}
					}
				else {
					PlotArray[j][i]=Average_Data[param[j]][i];
					PlotError[j][i]=Average_Error[param[j]][i];
					PlotCut[j][i]=Cut_Data[param[j]][i];
					}
				}
			}
		prev_xtype=type[0];
		prev_ytype=type[1];
		prev_xparam=param[0];
		prev_yparam=param[1];
		}		
	
	return 0;
}


int CVICALLBACK Add_All_Points (int panel, int control, int event,
		void *callbackData, int eventData1, int eventData2)
{
	int i;
	
	if(event==EVENT_COMMIT) {
		for(i=0;i<PlotLines;++i) {
			PlotCut[0][i]=0;
			PlotCut[1][i]=0;
			}
		fit_flag=FALSE;
		Plot_points(xplotstart,xplotstop,yplotstart,yplotstop);
		}
	return 0;
}


int CVICALLBACK CursorTimer (int panel, int control, int event,
		void *callbackData, int eventData1, int eventData2)
{   /* fix for displaying coordinates of points on initial "left click" 
	   rechecks cursor position for 5 seconds.
	*/
	double x,y;
	char graphlabel[100];
	
	if(event==EVENT_TIMER_TICK) {
		if(timercount>=TIMERSTOP) {
			SetCtrlAttribute(big,AP_CursorTimer,ATTR_ENABLED,0);
			return 0;
			}
		GetGraphCursor(big,AP_Grid,1,&x,&y);
		sprintf(graphlabel,"x = %g y = %g",x,y);
		SetCtrlAttribute(big, AP_Grid,ATTR_LABEL_TEXT,graphlabel);
		++timercount;
		}
	return 0;
}


int CVICALLBACK SavePlot (int panel, int control, int event,
		void *callbackData, int eventData1, int eventData2)
{
	int i,show;
	char filename[500],line[500];
	FILE *fp;
	
	if(event==EVENT_COMMIT) {
		GetCtrlVal(big,AP_plotname,filename);
		GetCtrlVal(big,AP_showcut,&show);
		fp=fopen(filename,"w");
		for(i=0;i<PlotLines;++i) {
			if((show && PlotCut[0][i]!=2 && PlotCut[1][i]!=2) ||
				(PlotCut[0][i]==0 && PlotCut[1][i]==0)) {
				sprintf(line,"%20.12g%20.12g%20.12g%20.12g\n",PlotArray[0][i],PlotArray[1][i],PlotError[0][i],PlotError[1][i]);
				fputs(line,fp);
				}
			}
		fclose(fp);
		}
		
	return 0;
}


int CVICALLBACK ZeroErr (int panel, int control, int event,
		void *callbackData, int eventData1, int eventData2)
{
	int i;
	
	if(event==EVENT_COMMIT) {
		for(i=0;i<2;++i)
			memset(PlotError[i],0,PlotLines*sizeof(double));
		}
	return 0;
}

int CVICALLBACK AutoGraphCtrl (int panel, int control, int event,
		void *callbackData, int eventData1, int eventData2)
{   // when autograph button is "on" the graph updates automatically when parameters are changed
	
	if(event==EVENT_COMMIT) 
	switch (control) {
		case AP_AUTOGRAPH:
			GetCtrlVal(big,AP_AUTOGRAPH,&autograph);
			break;
		case AP_AUTOFIT:
			GetCtrlVal(big,AP_AUTOFIT,&autofit);
			break;
		}
	return 0;
}

int CVICALLBACK PlotParam (int panel, int control, int event,
		void *callbackData, int eventData1, int eventData2)
{
	if(event==EVENT_VAL_CHANGED) {
		if(autograph) PlotSet(big,AP_Graph_button,EVENT_COMMIT,0,0,0);
		if(autofit) Draw_Fit(big,AP_LinearFit,EVENT_COMMIT,0,0,0);
		}
	return 0;
}
