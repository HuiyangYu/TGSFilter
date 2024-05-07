#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <functional>

using namespace std;

struct LenDisData
{
    vector<int> x;
    vector<uint64_t> y;
};

struct QualDisData
{
    vector<int> x;
    vector<float> y;
};


struct LineYData
{
    string bases;
    vector<float> data;
};

struct LinePlotData
{
    vector<int> x;
    vector<LineYData> y;
};

struct Selection
{
    string title;
    vector<string> ids;
};

struct OnePlotData
{
    LenDisData downLenDis;
    QualDisData downQualDis;
    LinePlotData downReadsQual;
    LinePlotData down5pReadsQual;
    LinePlotData down3pReadsQual;
    LinePlotData downBasesContents;
    LinePlotData down5pBasesContents;
    LinePlotData down3pBasesContents;
};

struct TwoPlotData
{
    LenDisData rawLenDis;
    LenDisData cleanLenDis;
    QualDisData rawQualDis;
    QualDisData cleanQualDis;
    LinePlotData rawReadsQual;
    LinePlotData cleanReadsQual;
    LinePlotData raw5pReadsQual;
    LinePlotData clean5pReadsQual;
    LinePlotData raw3pReadsQual;
    LinePlotData clean3pReadsQual;
    LinePlotData rawBasesContents;
    LinePlotData cleanBasesContents;
    LinePlotData raw5pBasesContents;
    LinePlotData clean5pBasesContents;
    LinePlotData raw3pBasesContents;
    LinePlotData clean3pBasesContents;
};

struct ThreePlotData
{
    LenDisData downLenDis;
    LenDisData cleanLenDis;
    LenDisData rawLenDis;
    QualDisData downQualDis;
    QualDisData rawQualDis;
    QualDisData cleanQualDis;
    LinePlotData downReadsQual;
    LinePlotData down5pReadsQual;
    LinePlotData down3pReadsQual;
    LinePlotData downBasesContents;
    LinePlotData down5pBasesContents;
    LinePlotData down3pBasesContents;
    LinePlotData rawReadsQual;
    LinePlotData cleanReadsQual;
    LinePlotData raw5pReadsQual;
    LinePlotData clean5pReadsQual;
    LinePlotData raw3pReadsQual;
    LinePlotData clean3pReadsQual;
    LinePlotData rawBasesContents;
    LinePlotData cleanBasesContents;
    LinePlotData raw5pBasesContents;
    LinePlotData clean5pBasesContents;
    LinePlotData raw3pBasesContents;
    LinePlotData clean3pBasesContents;
};
