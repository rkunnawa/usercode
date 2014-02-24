#include <iostream>
#include <iomanip>
#include <fstream>
#include <utility>
#include <TRandom.h>
#include <TH1D.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TROOT.h>
#include <TStyle.h>
#include "TCanvas.h"
#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TBox.h>
#include "TF1.h"
#include "TDirectory.h"
#include "TDirectoryFile.h"
#include "TGraph.h"
#include "TMath.h"
#include "TLegend.h"
#include "TLine.h"
#include "TLatex.h"

using namespace std;

// Truth binning

static const int nbins_truth = 1000;
static const double boundaries_truth[nbins_truth+1] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 180, 181, 182, 183, 184, 185, 186, 187, 188, 189, 190, 191, 192, 193, 194, 195, 196, 197, 198, 199, 200, 201, 202, 203, 204, 205, 206, 207, 208, 209, 210, 211, 212, 213, 214, 215, 216, 217, 218, 219, 220, 221, 222, 223, 224, 225, 226, 227, 228, 229, 230, 231, 232, 233, 234, 235, 236, 237, 238, 239, 240, 241, 242, 243, 244, 245, 246, 247, 248, 249, 250, 251, 252, 253, 254, 255, 256, 257, 258, 259, 260, 261, 262, 263, 264, 265, 266, 267, 268, 269, 270, 271, 272, 273, 274, 275, 276, 277, 278, 279, 280, 281, 282, 283, 284, 285, 286, 287, 288, 289, 290, 291, 292, 293, 294, 295, 296, 297, 298, 299, 300, 301, 302, 303, 304, 305, 306, 307, 308, 309, 310, 311, 312, 313, 314, 315, 316, 317, 318, 319, 320, 321, 322, 323, 324, 325, 326, 327, 328, 329, 330, 331, 332, 333, 334, 335, 336, 337, 338, 339, 340, 341, 342, 343, 344, 345, 346, 347, 348, 349, 350, 351, 352, 353, 354, 355, 356, 357, 358, 359, 360, 361, 362, 363, 364, 365, 366, 367, 368, 369, 370, 371, 372, 373, 374, 375, 376, 377, 378, 379, 380, 381, 382, 383, 384, 385, 386, 387, 388, 389, 390, 391, 392, 393, 394, 395, 396, 397, 398, 399, 400, 401, 402, 403, 404, 405, 406, 407, 408, 409, 410, 411, 412, 413, 414, 415, 416, 417, 418, 419, 420, 421, 422, 423, 424, 425, 426, 427, 428, 429, 430, 431, 432, 433, 434, 435, 436, 437, 438, 439, 440, 441, 442, 443, 444, 445, 446, 447, 448, 449, 450, 451, 452, 453, 454, 455, 456, 457, 458, 459, 460, 461, 462, 463, 464, 465, 466, 467, 468, 469, 470, 471, 472, 473, 474, 475, 476, 477, 478, 479, 480, 481, 482, 483, 484, 485, 486, 487, 488, 489, 490, 491, 492, 493, 494, 495, 496, 497, 498, 499, 500, 501, 502, 503, 504, 505, 506, 507, 508, 509, 510, 511, 512, 513, 514, 515, 516, 517, 518, 519, 520, 521, 522, 523, 524, 525, 526, 527, 528, 529, 530, 531, 532, 533, 534, 535, 536, 537, 538, 539, 540, 541, 542, 543, 544, 545, 546, 547, 548, 549, 550, 551, 552, 553, 554, 555, 556, 557, 558, 559, 560, 561, 562, 563, 564, 565, 566, 567, 568, 569, 570, 571, 572, 573, 574, 575, 576, 577, 578, 579, 580, 581, 582, 583, 584, 585, 586, 587, 588, 589, 590, 591, 592, 593, 594, 595, 596, 597, 598, 599, 600, 601, 602, 603, 604, 605, 606, 607, 608, 609, 610, 611, 612, 613, 614, 615, 616, 617, 618, 619, 620, 621, 622, 623, 624, 625, 626, 627, 628, 629, 630, 631, 632, 633, 634, 635, 636, 637, 638, 639, 640, 641, 642, 643, 644, 645, 646, 647, 648, 649, 650, 651, 652, 653, 654, 655, 656, 657, 658, 659, 660, 661, 662, 663, 664, 665, 666, 667, 668, 669, 670, 671, 672, 673, 674, 675, 676, 677, 678, 679, 680, 681, 682, 683, 684, 685, 686, 687, 688, 689, 690, 691, 692, 693, 694, 695, 696, 697, 698, 699, 700, 701, 702, 703, 704, 705, 706, 707, 708, 709, 710, 711, 712, 713, 714, 715, 716, 717, 718, 719, 720, 721, 722, 723, 724, 725, 726, 727, 728, 729, 730, 731, 732, 733, 734, 735, 736, 737, 738, 739, 740, 741, 742, 743, 744, 745, 746, 747, 748, 749, 750, 751, 752, 753, 754, 755, 756, 757, 758, 759, 760, 761, 762, 763, 764, 765, 766, 767, 768, 769, 770, 771, 772, 773, 774, 775, 776, 777, 778, 779, 780, 781, 782, 783, 784, 785, 786, 787, 788, 789, 790, 791, 792, 793, 794, 795, 796, 797, 798, 799, 800, 801, 802, 803, 804, 805, 806, 807, 808, 809, 810, 811, 812, 813, 814, 815, 816, 817, 818, 819, 820, 821, 822, 823, 824, 825, 826, 827, 828, 829, 830, 831, 832, 833, 834, 835, 836, 837, 838, 839, 840, 841, 842, 843, 844, 845, 846, 847, 848, 849, 850, 851, 852, 853, 854, 855, 856, 857, 858, 859, 860, 861, 862, 863, 864, 865, 866, 867, 868, 869, 870, 871, 872, 873, 874, 875, 876, 877, 878, 879, 880, 881, 882, 883, 884, 885, 886, 887, 888, 889, 890, 891, 892, 893, 894, 895, 896, 897, 898, 899, 900, 901, 902, 903, 904, 905, 906, 907, 908, 909, 910, 911, 912, 913, 914, 915, 916, 917, 918, 919, 920, 921, 922, 923, 924, 925, 926, 927, 928, 929, 930, 931, 932, 933, 934, 935, 936, 937, 938, 939, 940, 941, 942, 943, 944, 945, 946, 947, 948, 949, 950, 951, 952, 953, 954, 955, 956, 957, 958, 959, 960, 961, 962, 963, 964, 965, 966, 967, 968, 969, 970, 971, 972, 973, 974, 975, 976, 977, 978, 979, 980, 981, 982, 983, 984, 985, 986, 987, 988, 989, 990, 991, 992, 993, 994, 995, 996, 997, 998, 999, 1000
};

/*
static const int nbins_truth = 40;
static const double boundaries_truth[nbins_truth+1] = {
	100, 110, 120, 130, 140, 150, 160, 170,
	180, 190, 200, 210, 220, 230, 240, 250, 260, 270, 280, 290, 300, 310,
	320, 330, 340, 350, 360, 370, 380, 390, 400, 410, 420, 430, 440, 450,
	460, 470, 480, 490, 500
};
*/
// Measurement binning

static const int nbins_rec = 1000;
static const double boundaries_rec[nbins_rec+1] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 180, 181, 182, 183, 184, 185, 186, 187, 188, 189, 190, 191, 192, 193, 194, 195, 196, 197, 198, 199, 200, 201, 202, 203, 204, 205, 206, 207, 208, 209, 210, 211, 212, 213, 214, 215, 216, 217, 218, 219, 220, 221, 222, 223, 224, 225, 226, 227, 228, 229, 230, 231, 232, 233, 234, 235, 236, 237, 238, 239, 240, 241, 242, 243, 244, 245, 246, 247, 248, 249, 250, 251, 252, 253, 254, 255, 256, 257, 258, 259, 260, 261, 262, 263, 264, 265, 266, 267, 268, 269, 270, 271, 272, 273, 274, 275, 276, 277, 278, 279, 280, 281, 282, 283, 284, 285, 286, 287, 288, 289, 290, 291, 292, 293, 294, 295, 296, 297, 298, 299, 300, 301, 302, 303, 304, 305, 306, 307, 308, 309, 310, 311, 312, 313, 314, 315, 316, 317, 318, 319, 320, 321, 322, 323, 324, 325, 326, 327, 328, 329, 330, 331, 332, 333, 334, 335, 336, 337, 338, 339, 340, 341, 342, 343, 344, 345, 346, 347, 348, 349, 350, 351, 352, 353, 354, 355, 356, 357, 358, 359, 360, 361, 362, 363, 364, 365, 366, 367, 368, 369, 370, 371, 372, 373, 374, 375, 376, 377, 378, 379, 380, 381, 382, 383, 384, 385, 386, 387, 388, 389, 390, 391, 392, 393, 394, 395, 396, 397, 398, 399, 400, 401, 402, 403, 404, 405, 406, 407, 408, 409, 410, 411, 412, 413, 414, 415, 416, 417, 418, 419, 420, 421, 422, 423, 424, 425, 426, 427, 428, 429, 430, 431, 432, 433, 434, 435, 436, 437, 438, 439, 440, 441, 442, 443, 444, 445, 446, 447, 448, 449, 450, 451, 452, 453, 454, 455, 456, 457, 458, 459, 460, 461, 462, 463, 464, 465, 466, 467, 468, 469, 470, 471, 472, 473, 474, 475, 476, 477, 478, 479, 480, 481, 482, 483, 484, 485, 486, 487, 488, 489, 490, 491, 492, 493, 494, 495, 496, 497, 498, 499, 500, 501, 502, 503, 504, 505, 506, 507, 508, 509, 510, 511, 512, 513, 514, 515, 516, 517, 518, 519, 520, 521, 522, 523, 524, 525, 526, 527, 528, 529, 530, 531, 532, 533, 534, 535, 536, 537, 538, 539, 540, 541, 542, 543, 544, 545, 546, 547, 548, 549, 550, 551, 552, 553, 554, 555, 556, 557, 558, 559, 560, 561, 562, 563, 564, 565, 566, 567, 568, 569, 570, 571, 572, 573, 574, 575, 576, 577, 578, 579, 580, 581, 582, 583, 584, 585, 586, 587, 588, 589, 590, 591, 592, 593, 594, 595, 596, 597, 598, 599, 600, 601, 602, 603, 604, 605, 606, 607, 608, 609, 610, 611, 612, 613, 614, 615, 616, 617, 618, 619, 620, 621, 622, 623, 624, 625, 626, 627, 628, 629, 630, 631, 632, 633, 634, 635, 636, 637, 638, 639, 640, 641, 642, 643, 644, 645, 646, 647, 648, 649, 650, 651, 652, 653, 654, 655, 656, 657, 658, 659, 660, 661, 662, 663, 664, 665, 666, 667, 668, 669, 670, 671, 672, 673, 674, 675, 676, 677, 678, 679, 680, 681, 682, 683, 684, 685, 686, 687, 688, 689, 690, 691, 692, 693, 694, 695, 696, 697, 698, 699, 700, 701, 702, 703, 704, 705, 706, 707, 708, 709, 710, 711, 712, 713, 714, 715, 716, 717, 718, 719, 720, 721, 722, 723, 724, 725, 726, 727, 728, 729, 730, 731, 732, 733, 734, 735, 736, 737, 738, 739, 740, 741, 742, 743, 744, 745, 746, 747, 748, 749, 750, 751, 752, 753, 754, 755, 756, 757, 758, 759, 760, 761, 762, 763, 764, 765, 766, 767, 768, 769, 770, 771, 772, 773, 774, 775, 776, 777, 778, 779, 780, 781, 782, 783, 784, 785, 786, 787, 788, 789, 790, 791, 792, 793, 794, 795, 796, 797, 798, 799, 800, 801, 802, 803, 804, 805, 806, 807, 808, 809, 810, 811, 812, 813, 814, 815, 816, 817, 818, 819, 820, 821, 822, 823, 824, 825, 826, 827, 828, 829, 830, 831, 832, 833, 834, 835, 836, 837, 838, 839, 840, 841, 842, 843, 844, 845, 846, 847, 848, 849, 850, 851, 852, 853, 854, 855, 856, 857, 858, 859, 860, 861, 862, 863, 864, 865, 866, 867, 868, 869, 870, 871, 872, 873, 874, 875, 876, 877, 878, 879, 880, 881, 882, 883, 884, 885, 886, 887, 888, 889, 890, 891, 892, 893, 894, 895, 896, 897, 898, 899, 900, 901, 902, 903, 904, 905, 906, 907, 908, 909, 910, 911, 912, 913, 914, 915, 916, 917, 918, 919, 920, 921, 922, 923, 924, 925, 926, 927, 928, 929, 930, 931, 932, 933, 934, 935, 936, 937, 938, 939, 940, 941, 942, 943, 944, 945, 946, 947, 948, 949, 950, 951, 952, 953, 954, 955, 956, 957, 958, 959, 960, 961, 962, 963, 964, 965, 966, 967, 968, 969, 970, 971, 972, 973, 974, 975, 976, 977, 978, 979, 980, 981, 982, 983, 984, 985, 986, 987, 988, 989, 990, 991, 992, 993, 994, 995, 996, 997, 998, 999, 1000
	};

/*
 static const int nbins_rec = 40;
 static const double boundaries_rec[nbins_rec+1] = {
 100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200, 210, 220,
 230, 240, 250, 260, 270, 280, 290, 300, 310, 320, 330, 340, 350, 360,
 370, 380, 390, 400, 410, 420, 430, 440, 450, 460, 470, 480, 490, 500
 };
*/

/*
static const int nbins_rec = 20;
static const double boundaries_rec[nbins_rec+1] = {
  0,5,10,15,20,25,30,35,40,45,50,55,60,70,80,90,110,130,170,250,300
};

static const int nbins_truth = 20;
static const double boundaries_truth[nbins_truth+1] = {
  0,5,10,15,20,25,30,35,40,45,50,55,60,70,80,90,110,130,170,250,300
};
*/
//changed it for yaxian's comparison

//static const int nbins_recrebin = 20;
//static const double boundaries_recrebin[nbins_recrebin+1] = {
//  0,5,10,15,20,25,30,35,40,45,50,55,60,70,80,90,110,130,170,250,330
//};

//temporary nbins_recrebin to accomodate the old binning style for the NLO comparison
static const int nbins_recrebin = 26;
static const double boundaries_recrebin[nbins_recrebin+1] = {
  30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,200,220,240,260,280,300,340,380,420,460,500
};


static const int nbins_recrebinM = 19;
static const double boundaries_recrebinM[nbins_recrebinM+1] = {
  0,5,10,15,20,25,30,35,40,45,50,55,60,70,80,90,110,130,170,250
};
//the original bin values. 
//static const int nbins_recrebin = 18;
//static const double boundaries_recrebin[nbins_recrebin+1] = {
//	30, 40, 50, 60, 70, 80, 90,100, 110, 120, 130, 140, 150, 160,
//	170, 180, 200, 240, 300
//};

//static const int nbins_recrebinM = 18;
//static const double boundaries_recrebinM[nbins_recrebinM+1] = {
//	30, 40, 50, 60, 70, 80, 90,100, 110, 120, 130, 140, 150, 160,
//	170, 180, 200, 240, 300
//};

static const int nbins_recrebin_Npart = 1;
static const double boundaries_recrebin_Npart[nbins_recrebin_Npart+1] = {
  100, 300//changed from 100 to 30
};

static const int nColor = 5;
static const int colorCode[nColor] = {
	1, 2, kGreen+1, 4, 6
};

// Algos, following Ying's convension
static const int nAlgos = 8;
static const int BinLabelN = 11;
static const char *algoName[nAlgos] = { "", "icPu5", "akPu2PF", "ak3PF", "akPu4PF", "akPu2Calo", "akPu3Calo", "akPu4Calo" };
static const char *algoNamePP[nAlgos] = { "", "icPu5", "ak2PF", "ak3PF", "ak4PF", "ak2Calo", "ak3Calo", "ak4Calo" };
static const char *algoNameGen[nAlgos] = { "", "icPu5", "akPu2PF", "ak3PF", "akPu4PF", "akPu2PF", "akPu3PF", "akPu4PF" };
static const char *BinLabel[BinLabelN] = {"100-110", "110-120", "120-130", "130-140", "140-150", "150-160", "160-170", "170-180", "180-200", "200-240","240-300" };

// Centrality binning
//const int nbins_cent=	     6;
//Double_t boundaries_cent[nbins_cent+1] = {   0,2,4,12,20,28,36   };// multiply by 2.5 to get your actual centrality % (old 2011 data) 
//Double_t ncoll[nbins_cent] = { 1660, 1310, 745, 251, 62.8, 10.8 };
//the following is due to request form Yaxian 
const int nbins_cent = 1;
Double_t boundaries_cent[nbins_cent+1] = {0,40};
Double_t ncoll[nbins_cent] = {362.24}; //use taa instead of ncoll. 
//Double_t TAA[nbins_cent] = {};

// TAA 

TGraphErrors *tTAAerr[6]={0};
TGraphErrors *tTAAerrNpart=0;


//static const double trigEffInc[nbins_recrebin] = {
//  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.966051, 0.99591, 0.999665, 1., 1., 1., 1.
//};


// make a histogram from TF1 function
TH1F *functionHist(TF1 *f, TH1F* h,char *fHistname)
{
	TH1F *hF = (TH1F*)h->Clone(fHistname);
	for (int i=1;i<=h->GetNbinsX();i++)
	{
		double var = f->Integral(h->GetBinLowEdge(i),h->GetBinLowEdge(i+1))/h->GetBinWidth(i);
		hF->SetBinContent(i,var);
		hF->SetBinError(i,0);
	}
	return hF;
}

TLegend *myLegend(double x1,double y1,double x2, double y2)
{
	TLegend *leg = new TLegend(x1,y1,x2,y2);
	leg->SetBorderSize(0);
	leg->SetFillStyle(0);
	return leg; 
	
}


// draw envelope using a systematic uncertainty histogram
TH1F* drawEnvelope(TH1F *h,char *opt,int color = kGray,int fillStyle = 0, int fillColor = 0,double shift = 0)
{
	TH1F *hClone = (TH1F*) h->Clone(Form("%s_mirror",h->GetTitle()));
	TH1F *hMirror = (TH1F*) h->Clone(Form("%s_mirror2",h->GetTitle()));
	for (int i=1;i<=h->GetNbinsX();i++)
	{
		double val = h->GetBinContent(i);
		hMirror->SetBinContent(i,1-fabs(val-1)+shift);
		hClone->SetBinContent(i,val+shift);
	}
	
	//   hMirror->SetLineStyle(2);
	//   h->SetLineStyle(2);
	hMirror->SetLineColor(color);
	hMirror->SetFillColor(fillColor);
	hMirror->SetFillStyle(fillStyle);
	hClone->SetLineColor(color);
	hClone->SetFillColor(fillColor);
	hClone->SetFillStyle(fillStyle);
	hClone->Draw(opt);
	hMirror->Draw(opt);
	return hMirror;
}

void makeHistTitle(TH1 *h,char *title, char *xTitle, char *yTitle, int color = -1, bool centerTitle = 1)
{
	h->SetTitle(title);
	h->SetXTitle(xTitle);
	h->SetYTitle(yTitle);
	
	if (centerTitle) {
		h->GetXaxis()->CenterTitle();
		h->GetYaxis()->CenterTitle();
		
	}
	
	if (color!=-1) {
		h->SetLineColor(color);
		h->SetMarkerColor(color);
	}
	
	h->GetYaxis()->SetNdivisions(610); 
	h->GetXaxis()->SetNdivisions(505);

	
	h->GetYaxis()->SetLabelFont(43);
	h->GetYaxis()->SetTitleFont(43);
	h->GetYaxis()->SetLabelSize(20);
	h->GetYaxis()->SetTitleSize(22);
	h->GetYaxis()->SetTitleOffset(2.6);
	
	h->GetXaxis()->SetLabelFont(43);
	h->GetXaxis()->SetTitleFont(43);
	h->GetXaxis()->SetLabelSize(20);
	h->GetXaxis()->SetTitleSize(22);
	h->GetXaxis()->SetTitleOffset(3.1);
	
	h->GetXaxis()->SetNoExponent();
	h->GetXaxis()->SetMoreLogLabels();
	
	h->GetXaxis()->SetTitleOffset(2.4);
	
}



class SysData
{
public:
	SysData() {
		for (int i=0;i<=nbins_cent;i++) {
			hSys[i]     = new TH1F(Form("hSys_cent%d",i), Form("Totalsys_cent%d",i),nbins_recrebin, boundaries_recrebin);
			hSysGeneral[i]= new TH1F(Form("hSysGeneral_cent%d",i), Form("TotalsysGeneral_cent%d",i),nbins_recrebin, boundaries_recrebin);
			hSysJEC[i]  = new TH1F(Form("hSysJEC_cent%d",i), Form("JECsys_cent%d",i),nbins_recrebin, boundaries_recrebin);
			hSysEff[i]  = new TH1F(Form("hSysEff_cent%d",i), Form("Effsys_cent%d",i),nbins_recrebin, boundaries_recrebin);
			hSysSmear[i]  = new TH1F(Form("hSysSmear_cent%d",i), Form("Smearsys_cent%d",i),nbins_recrebin, boundaries_recrebin);
			hSysIter[i] = new TH1F(Form("hSysIter_cent%d",i), Form("Itersys_cent%d",i),nbins_recrebin, boundaries_recrebin);
			hSys[i]->SetLineColor(kGray);
			hSysJEC[i]->SetLineColor(4);
			hSysSmear[i]->SetLineColor(kGreen+1);
			hSysIter[i]->SetLineColor(2);
		}  
	}
	TH1F *hSys[nbins_cent+1];
	TH1F *hSysGeneral[nbins_cent+1];
	TH1F *hSysJEC[nbins_cent+1];
	TH1F *hSysEff[nbins_cent+1];
	TH1F *hSysSmear[nbins_cent+1];
	TH1F *hSysIter[nbins_cent+1];
	TH1F *hSysNoise[nbins_cent+1];
	
	void calcTotalSys(int i) {
		TF1 *fNoise = new TF1("f","1+0.3*0.16*abs(1-([0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x))");
		fNoise->SetParameters(0.9521,0.001105,-9.397e-6,3.32e-8,-5.618e-11);
		hSysNoise[i] = functionHist(fNoise,hSys[i],Form("hSysNoise_cent%d",i));
		hSysNoise[i]->SetName(Form("hSysNoise_cent%d",i));
		hSysNoise[i]->SetLineColor(6);
		for (int j=1;j<=hSys[i]->GetNbinsX();j++) {
			double effSys = 0.01;
			hSysEff[i]->SetBinContent(j,1+effSys);
			hSysSmear[i]->SetBinContent(j,1.02);
			double JECSys = hSysJEC[i]->GetBinContent(j)-1;
			double SmearSys = hSysSmear[i]->GetBinContent(j)-1;
			double IterSys = hSysIter[i]->GetBinContent(j)-1; 
			double NoiseSys = hSysNoise[i]->GetBinContent(j)-1; 
			cout <<effSys<<" "<<JECSys<<" "<<IterSys<<endl;
			double totalSys = sqrt( effSys * effSys +
								   JECSys * JECSys +
								   SmearSys * SmearSys +
								   NoiseSys * NoiseSys +
								   IterSys* IterSys
								   );
			hSys[i]->SetBinContent(j,totalSys+1);	
			hSys[i]->SetLineWidth(2);				
		}
	}
	
	void calcTotalSysNoUnfolding(int i) {
		TF1 *fNoise = new TF1("f","1+0.3*0.16*abs(1-([0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x))");
		fNoise->SetParameters(0.9521,0.001105,-9.397e-6,3.32e-8,-5.618e-11);
		hSysNoise[i] = functionHist(fNoise,hSys[i],Form("hSysNoise_cent%d",i));
		hSysNoise[i]->SetName(Form("hSysNoise_cent%d",i));
		hSysNoise[i]->SetLineColor(6);
		for (int j=1;j<=hSysGeneral[i]->GetNbinsX();j++) {
			double effSys = 0.01;
			hSysEff[i]->SetBinContent(j,1+effSys);
			hSysSmear[i]->SetBinContent(j,1.02);
			double JECSys = hSysJEC[i]->GetBinContent(j)-1;
			double SmearSys = hSysSmear[i]->GetBinContent(j)-1;
			double NoiseSys = hSysNoise[i]->GetBinContent(j)-1; 
			double totalSys = sqrt( effSys * effSys +
								   JECSys * JECSys +
								   SmearSys * SmearSys +
								   NoiseSys * NoiseSys
								   );
			hSysGeneral[i]->SetBinContent(j,totalSys+1);	
			hSysGeneral[i]->SetLineWidth(2);				
		}
	}
	
	void Draw(TH1F *h,int i) {
		for (int j=1;j<=hSys[i]->GetNbinsX();j++) {
			double val = h->GetBinContent(j);
			double err = hSys[i]->GetBinContent(j)-1;
			cout << "Sys Value Check" <<val<<" "<<err<<" "<<h->GetBinLowEdge(j)<<" "<<val*(1-err)<<" "<<h->GetBinLowEdge(j+1)<<" "<<val*(1+err)<<endl;
			TBox *b = new TBox(h->GetBinLowEdge(j),val*(1-err),h->GetBinLowEdge(j+1),val*(1+err));
			//b->SetFillColor(kGray);
			b->SetFillStyle(1001);
			//b->SetLineColor(kGray);
			
			
			//***********For Gunther's Color Systematics Band Peference
			b->SetFillColor(5);
			b->SetLineColor(5);
			b->Draw();
		}
	}
	
	void DrawTGraph(TGraphErrors *h,int i) {
			double xv;
			double val;
			for(int j=0;j<h->GetN();j++){	
			h->GetPoint(j,xv,val);
			double err = hSysGeneral[i]->GetBinContent(j+1)-1;
			//cout <<"value" <<val<<" "<<err<<" "<<hSysGeneral[i]->GetBinLowEdge(j+1)<<" "<<val*(1-err)<<" "<<hSysGeneral[i]->GetBinLowEdge(j+2)<<" "<<val*(1+err)<<endl;
			TBox *b_ = new TBox(hSysGeneral[i]->GetBinLowEdge(j+1),val*(1-err),hSysGeneral[i]->GetBinLowEdge(j+2),val*(1+err));
			b_->SetFillColor(kGray);
			b_->SetFillStyle(1001);
			b_->SetLineColor(kGray);
			b_->Draw();
		}
	}
	
		void DrawUnfoErr(TH1F *h,int i) {
			for (int j=1;j<=hSysIter[i]->GetNbinsX();j++) {
			double val = h->GetBinContent(j);
			double err = hSysIter[i]->GetBinContent(j)-1;
			//cout <<"value" << val<<" "<<err<<" "<<h->GetBinLowEdge(j)<<" "<<val*(1-err)<<" "<<h->GetBinLowEdge(j+1)<<" "<<val*(1+err)<<endl;
			TBox *b = new TBox(h->GetBinLowEdge(j),val*(1-err),h->GetBinLowEdge(j+1),val*(1+err));
			b->SetFillColor(29);
			b->SetFillStyle(1001);
			b->SetLineColor(29);
			b->Draw();
		}
	}
	
	
	void DrawNpartSys(double yvalue,int i,double xvalue) {
			double yerrorNpart[6]= {0.0409, 0.0459,0.0578,0.0944, 0.143, 0.176 };
			double err = hSys[i]->GetBinContent(1)-1;
			TBox *b = new TBox(xvalue-6.,yvalue*(1-err-yerrorNpart[i]),xvalue+6.,yvalue*(1+err+yerrorNpart[i]));
			//cout << "value " << yvalue<<" err   "<<err<<" xvalue  "<<xvalue<<" "<<yvalue*(1-err)<<" "<<yvalue*(1+err)<<endl;
			b->SetFillColor(kGray);
			b->SetFillStyle(1001);
			b->SetLineColor(kGray);
			
			//***********For Gunther's Color Systematics Band Peference
			//b->SetFillColor(5);
			//b->SetLineColor(5);
			
			b->Draw();
		
	}
	
	
	void DrawComponent(int i) {
		calcTotalSys(i);
		TH1D *h = new TH1D(Form("hSysTmp_cent%d",i),"",nbins_recrebin, boundaries_recrebin);
		makeHistTitle(h,"","Jet p_{T} (GeV/c)","Systematic uncertainty");
		h->SetAxisRange(-0.25,0.4,"Y");
		h->Draw();
		TH1F* sys = drawEnvelope(hSys[i],"same",hSys[i]->GetLineColor(),1001,hSys[i]->GetLineColor(),-1);
		TH1F* sysIter = drawEnvelope(hSysIter[i],"same",hSysIter[i]->GetLineColor(),3004,hSysIter[i]->GetLineColor(),-1);
		TH1F* sysJEC = drawEnvelope(hSysJEC[i],"same",hSysJEC[i]->GetLineColor(),3005,hSysJEC[i]->GetLineColor(),-1);
		TH1F* sysSmear =  drawEnvelope(hSysSmear[i],"same",hSysSmear[i]->GetLineColor(),3001,hSysSmear[i]->GetLineColor(),-1);
		TH1F* sysEff = drawEnvelope(hSysEff[i],"same",hSysEff[i]->GetLineColor(),3002,hSysEff[i]->GetLineColor(),-1);
		TH1F* sysNoise = drawEnvelope(hSysNoise[i],"same",hSysNoise[i]->GetLineColor(),3001,hSysNoise[i]->GetLineColor(),-1);
		TLine *l = new TLine(h->GetBinLowEdge(1),0,h->GetBinLowEdge(h->GetNbinsX()+1),0);
		l->Draw();
		TLine *l2 = new TLine(h->GetBinLowEdge(1),-0.25,h->GetBinLowEdge(1),0.4);
		l2->Draw();
		TLegend *leg = myLegend(0.52,0.6,0.95,0.93);
		leg->SetTextSize(0.043);
		leg->AddEntry(sys,"Total Systematics","f");
		leg->AddEntry(sysIter,"Unfolding","f");
		leg->AddEntry(sysJEC,"Jet Energy Scale","f");
		leg->AddEntry(sysEff,"Jet Trigger Efficiency","f");
		leg->AddEntry(sysSmear,"UE fluctuation","f");
		leg->AddEntry(sysNoise,"HCAL Noise","f");
		if (i==nbins_cent-1)leg->Draw();
	}
	
 
};

class JetData
{
public:
	JetData(char *fileName, char *jetTree, char *genJetTree, bool loadGenJet = 1) {
		cout <<"Open "<<fileName<<endl;
		tFile = new TFile(fileName,"read");
		tEvt = (TTree*)tFile->Get("hiEvtAnalyzer/HiTree");
		tJet = (TTree*)tFile->Get(jetTree);
		tJet->SetBranchAddress("jtpt" , jtpt );
		tJet->SetBranchAddress("trackMax" , trackMax );
		tJet->SetBranchAddress("chargedMax",chargedMax);
		tJet->SetBranchAddress("refpt", refpt);
		tJet->SetBranchAddress("nref" ,&njets);
		tJet->SetBranchAddress("jteta", jteta);
		tJet->SetBranchAddress("pthat",&pthat);
		if (loadGenJet) tGenJet = (TTree*)tFile->Get(genJetTree);
		if (loadGenJet) tGenJet->SetBranchAddress("ngen" ,&ngen);
		if (loadGenJet) tGenJet->SetBranchAddress("genpt", genpt);
		if (loadGenJet) tGenJet->SetBranchAddress("gensubid", gensubid);
		tEvt->SetBranchAddress("hiBin",&bin  );
		tEvt->SetBranchAddress("vz",&vz  );
		tJet->AddFriend(tEvt);
	};
	TFile *tFile;
	TTree *tJet;
	TTree *tGenJet;
	TTree *tEvt;
	float jtpt[1000];
	float refpt[1000];
	float jteta[1000];
	float trackMax[1000];
	float chargedMax[1000];
	float genpt[1000];
	int gensubid[1000];
	float vz;
	float pthat;
	int njets;
	int ngen;
	int bin;      
};

class UnfoldingHistos
{
public:
	UnfoldingHistos(int i) {
		hResTrue       = new TH1F(Form("Restrue_cent%d",i), Form(" RooUnfold Reconstructed Truth cent%d",i)    ,    nbins_truth, boundaries_truth); //Gen
		hResMeas       = new TH1F(Form("Resmeas_cent%d",i), Form(" RooUnfold ResMeasured_cent1%d",i),    nbins_rec, boundaries_rec); //Reco	
		hGen           = new TH1F(Form("hGen_cent%d",i) , Form(" Generator Truth_cent%d",i)    ,    nbins_truth, boundaries_truth); //Gen
		hRecoMC        = new TH1F(Form("hRecoMC_cent%d",i) , Form(" Generator Reco_cent%d",i)    ,    nbins_truth, boundaries_truth); //RecoMC
		hUnfoMC        = new TH1F(Form("hUnfoMC_cent%d",i),Form("Generator Reco Unfolded cent%d",i), nbins_truth, boundaries_truth); //unfolded reco MC
		hMatrix        = new TH2F(Form("hMatrix_cent%d",i), Form(" Response Matrix;Genjet^{p_{T}};Recojet^{p_{T}}",i), nbins_truth, boundaries_truth,nbins_rec, boundaries_rec); //Matrix
		hMatrixRebin   = new TH2F(Form("hMatrixRebin_cent%d",i), Form(" Response Matrix;Genjet^{p_{T}};Recojet^{p_{T}}",i), nbins_recrebinM, boundaries_recrebinM,nbins_recrebinM, boundaries_recrebinM); //Matrix Rebin
		hMeas          = new TH1F(Form("hMeas_cent%d",i)       , Form("Measured jet p_{T} spectra_cent%d",i)       ,nbins_rec,boundaries_rec);
		//hCombined      = new TH1F(Form("hCombined_cent%d",i)   , Form("Combined measured jet p_{T} spectra_cent%d",i)       ,nbins_rec,boundaries_rec);
		//hMeas_80        = new TH1F(Form("hMeas_80_cent%d",i)       , Form("Measured jet p_{T} from HLT80 spectra_cent%d",i)       ,nbins_rec,boundaries_rec);
		//hMeas_65        = new TH1F(Form("hMeas_65_cent%d",i)       , Form("Measured jet p_{T} from HLT65 spectra_cent%d",i)       ,nbins_rec,boundaries_rec);
		//hMeas_55        = new TH1F(Form("hMeas_55_cent%d",i)       , Form("Measured jet p_{T} from HLT55 spectra_cent%d",i)       ,nbins_rec,boundaries_rec);
		//hMeas80JECSys    = new TH1F(Form("hMeas80JECSys_cent%d",i) , Form(" Measured jet p_{T} 80 sectra for JEC systematics_cent%d",i)       ,nbins_rec,boundaries_rec);
		//hCombinedJECSys    = new TH1F(Form("hCombinedJECSys_cent%d",i) , Form(" Combined Measured jet p_{T} sectra for JEC systematics_cent%d",i)       ,nbins_rec,boundaries_rec);
		//hMeas80SmearSys    = new TH1F(Form("hMeas80SmearSys_cent%d",i) , Form(" Measured jet p_{T}80 sectra for Smear systematics_cent%d",i)       ,nbins_rec,boundaries_rec);
		//hCombinedSmearSys    = new TH1F(Form("hCombinedSmearSys_cent%d",i) , Form(" Combined Measured jet p_{T} sectra for Smear systematics_cent%d",i)       ,nbins_rec,boundaries_rec);
		//hMeas65JECSys    = new TH1F(Form("hMeas65JECSys_cent%d",i) , Form(" Measured jet p_{T}65 sectra for JEC systematics_cent%d",i)       ,nbins_rec,boundaries_rec);
		//hMeas65SmearSys    = new TH1F(Form("hMeas65SmearSys_cent%d",i) , Form(" Measured jet p_{T}65 sectra for Smear systematics_cent%d",i)       ,nbins_rec,boundaries_rec);
		//hMeas55JECSys    = new TH1F(Form("hMeas55JECSys_cent%d",i) , Form(" Measured jet p_{T}55 sectra for JEC systematics_cent%d",i)       ,nbins_rec,boundaries_rec);
		//hMeas55SmearSys    = new TH1F(Form("hMeas55SmearSys_cent%d",i) , Form(" Measured jet p_{T}55 sectra for Smear systematics_cent%d",i)       ,nbins_rec,boundaries_rec);
		hMeasJECSys    = new TH1F(Form("hMeasJECSys_cent%d",i) , Form(" Measured jet p_{T} sectra for JEC systematics_cent%d",i)       ,nbins_rec,boundaries_rec);
		hMeasSmearSys    = new TH1F(Form("hMeasSmearSys_cent%d",i) , Form(" Measured jet p_{T} sectra for Smear systematics_cent%d",i)       ,nbins_rec,boundaries_rec);
		hMeas          ->Sumw2();
		/*
		hMeas_80        ->Sumw2();
		hMeas_65        ->Sumw2();
		hMeas_55        ->Sumw2();
		hCombined      ->Sumw2();
		hCombinedJECSys->Sumw2();
		hCombinedSmearSys->Sumw2();
		*/
		hMeasJECSys    ->Sumw2();
		hMeasSmearSys  ->Sumw2();
		/*
		hMeas80JECSys    ->Sumw2();
		hMeas80SmearSys  ->Sumw2();
		hMeas65JECSys    ->Sumw2();
		hMeas65SmearSys  ->Sumw2();
		hMeas55JECSys    ->Sumw2();
		hMeas55SmearSys  ->Sumw2();
		*/
		hResTrue       ->Sumw2();
		hGen           ->Sumw2();
		hRecoMC        ->Sumw2();
		hUnfoMC ->Sumw2();
		hResMeas       ->Sumw2();
		hMatrix        ->Sumw2();
		hMatrixRebin   ->Sumw2();
	};
	TH1F* hGen;
	//start adding by raghav
	TH1F* hRecoMC;
	TH1F* hUnfoMC;
	/*
	TH1F* hTrue;//to check with the pp data distribution
	TH1F* hUnfold;// make sure that unfolding works properly
	TH1F* hMC;
	//Try to define the different histograms to take in the entries for different HLT datasets here
	/*
	TH1F* hMeas_80;
	TH1F* hMeas_65;
	TH1F* hMeas_55;
	TH1F* hCombined;
        TH1F* hCombinedJECSys;
	TH1F* hCombinedSmearSys;
	TH1F* hMeas80JECSys;
	TH1F* hMeas80SmearSys;
	TH1F* hMeas65JECSys;
	TH1F* hMeas65SmearSys;
	TH1F* hMeas55JECSys;
	TH1F* hMeas55SmearSys;
	*/
	//end adding by raghav
	TH1F* hMCMeas;
	TH1F* hMeas;
	TH1F* hMeasJECSys;
	TH1F* hMeasSmearSys;
	TH1F* hResTrue;  
	TH1F* hResMeas;
	TH1F* hReco;
	TH1F* hRecoBinByBin;
	TH1F* hRecoJECSys;
	TH1F* hRecoSmearSys;
	TH1F* hRecoIterSys[10];  
	TH2F* hMatrix;
	TH2F* hMatrixRebin;
	TH2F* hMatrixNorm;
	TH2F* hResponse;
	TH2F* hResponseNorm;
	TH1F *hMeasMatch;
};

// Remove bins with error > central value
void cleanup(TH1F *h)
{
	for (int i=1;i<=h->GetNbinsX();i++)
	{
		double val1 = h->GetBinContent(i);
		double valErr1 = h->GetBinError(i);
		if (valErr1>=val1) {
			h->SetBinContent(i,0);
			h->SetBinError(i,0);
		}
	}   
	
}

// Remove error 
void removeError(TH1F *h)
{
	for (int i=1;i<=h->GetNbinsX();i++)
	{
		h->SetBinError(i,0);
	}   
	
}

// Remove Zero
void removeZero(TH1 *h)
{
  double min = 0;
  for(int i = 1;i<h->GetNbinsX();i++){
    if(h->GetBinContent(i)>min&&h->GetBinContent(i)>0)
      min = h->GetBinContent(i);
  }

  for(int i = 1;i<h->GetNbinsX();i++){
    if(h->GetBinContent(i) == 0){
      h->SetBinContent(i,min/10.);
      h->SetBinError(i,min/10.);
    }
  }
}


// rebin the spectra
TH1F *rebin(TH1F *h, char *histName)
{
	TH1F *hRebin = new TH1F(Form("%s_rebin",h->GetName()),Form("rebin %s",h->GetTitle()),nbins_recrebin,boundaries_recrebin);
	for (int i=1;i<=h->GetNbinsX();i++)
	{
		double val=h->GetBinContent(i);
		double valErr=h->GetBinError(i);
		int binNum = hRebin->FindBin(h->GetBinCenter(i));
		double val1 = hRebin->GetBinContent(binNum);
		double valErr1 = hRebin->GetBinError(binNum);
		hRebin->SetBinContent(binNum,val+val1);
		hRebin->SetBinError(binNum,sqrt(valErr1*valErr1+valErr*valErr));
	}
	cleanup(hRebin);
	hRebin->SetName(histName);
	return hRebin;
}

// rebin the spectra
TH1F *rebin2(TH1F *h, char *histName)
{
  TH1F *hRebin = new TH1F(Form("%s_rebin",h->GetName()),Form("rebin %s",h->GetTitle()),nbins_rec,boundaries_rec);
  for (int i=1;i<=h->GetNbinsX();i++)
    {
      double val=h->GetBinContent(i);
      double valErr=h->GetBinError(i);
      int binNum = hRebin->FindBin(h->GetBinCenter(i));
      double val1 = hRebin->GetBinContent(binNum);
      double valErr1 = hRebin->GetBinError(binNum);
      hRebin->SetBinContent(binNum,val+val1);
      hRebin->SetBinError(binNum,sqrt(valErr1*valErr1+valErr*valErr));
    }
  cleanup(hRebin);
  hRebin->SetName(histName);
  return hRebin;
}


// Npart rebin the spectra  
TH1F *rebin_Npart(TH1F *h, char *histName)
{
	TH1F *hRebin = new TH1F(Form("%s_rebin",h->GetName()),Form("rebin %s",h->GetTitle()),nbins_recrebin_Npart,boundaries_recrebin_Npart);
	for (int i=1;i<=h->GetNbinsX();i++)
	{
		double val=h->GetBinContent(i);
		double valErr=h->GetBinError(i);
		int binNum = hRebin->FindBin(h->GetBinCenter(i));
		double val1 = hRebin->GetBinContent(binNum);
		double valErr1 = hRebin->GetBinError(binNum);
		hRebin->SetBinContent(binNum,val+val1);
		hRebin->SetBinError(binNum,sqrt(valErr1*valErr1+valErr*valErr));
	}
	cleanup(hRebin);
	hRebin->SetName(histName);
	return hRebin;
}

// divide by bin width
void divideBinWidth(TH1 *h)
{
	h->Sumw2();
	for (int i=0;i<=h->GetNbinsX();i++)
	{
		Float_t val = h->GetBinContent(i);
		Float_t valErr = h->GetBinError(i);
		val/=h->GetBinWidth(i);
		valErr/=h->GetBinWidth(i);
		h->SetBinContent(i,val);
		h->SetBinError(i,valErr);
	}
	h->GetXaxis()->CenterTitle();
	h->GetYaxis()->CenterTitle();
}


// make systematic histogram
	 void checkMaximumSys(TH1F *hSys, TH1F *h, int opt=0,double minVal = 1)
{
	if (h->GetNbinsX()!=hSys->GetNbinsX()) {
		cout <<"ERROR! Different NBins in subroutine checkMaximumSys!"<<endl;
	} else {
		double val = minVal;
		for (int i=1;i<=h->GetNbinsX();i++) {
			//cout <<i<<" "<<val<<" "<<hSys->GetBinContent(i)<<" "<<h->GetBinContent(i)<<endl;
			if (h->GetBinContent(i)==0) continue;
			if (opt==0) val=minVal;
			if (fabs(hSys->GetBinContent(i))>val) val = fabs(hSys->GetBinContent(i));
			if (fabs(h->GetBinContent(i)-1)+1>val) val=fabs(h->GetBinContent(i)-1)+1;
			hSys->SetBinContent(i,val);
		}
	}
}



void makeMultiPanelCanvasWithGap(TCanvas*& canv,
								 const Int_t columns,
								 const Int_t rows,
								 const Float_t leftOffset,
								 const Float_t bottomOffset,
								 const Float_t leftMargin,
								 const Float_t bottomMargin,
								 const Float_t edge, const Float_t asyoffset) {
	if (canv==0) {
		Error("makeMultiPanelCanvas","Got null canvas.");
		return;
	}
	canv->Clear();
	
	TPad* pad[columns][rows];
	
	Float_t Xlow[columns];
	Float_t Xup[columns];
	Float_t Ylow[rows];
	Float_t Yup[rows];
	
	Float_t PadWidth =
	(1.0-leftOffset)/((1.0/(1.0-leftMargin)) +
					  (1.0/(1.0-edge))+(Float_t)columns-2.0);
	Float_t PadHeight =
	(1.0-bottomOffset)/((1.0/(1.0-bottomMargin)) +
						(1.0/(1.0-edge))+(Float_t)rows-2.0);
	
	//PadHeight = 0.5*PadHeight;
	
	Xlow[0] = leftOffset;
	Xup[0] = leftOffset + PadWidth/(1.0-leftMargin);
	Xup[columns-1] = 1;
	Xlow[columns-1] = 1.0-PadWidth/(1.0-edge);
	
	Yup[0] = 1;
	Ylow[0] = 1.0-PadHeight/(1.0-edge);
	Ylow[rows-1] = bottomOffset;
	Yup[rows-1] = bottomOffset + PadHeight/(1.0-bottomMargin);
	
	for(Int_t i=1;i<columns-1;i++) {
		Xlow[i] = Xup[0] + (i-1)*PadWidth;
		Xup[i] = Xup[0] + (i)*PadWidth;
	}
	Int_t ct = 0;
	for(Int_t i=rows-2;i>0;i--) {
		if(i==rows-2){
			Ylow[i] = Yup[rows-1] + ct*PadHeight;
			Yup[i] = Yup[rows-1] + (ct+1)*PadHeight;
		}else{
			Ylow[i] = Yup[rows-1] + ct*PadHeight;
			Yup[i] = Yup[rows-1] + (ct+1)*PadHeight;
			//Yup[i] = 0.2*Yup[i];
		}
		ct++;
	}
	
	TString padName;
	for(Int_t i=0;i<columns;i++) {
		for(Int_t j=0;j<rows;j++) {
			canv->cd();
			padName = Form("p_%d_%d",i,j);
			//pad[i][j] = new TPad(padName.Data(),padName.Data(),
			//Xlow[i],Ylow[j],Xup[i],Yup[j]);
			// this is hacked version to create aysmmetric pads around low 
			if(j==0){
				pad[i][j] = new TPad(padName.Data(),padName.Data(),
									 Xlow[i],Ylow[j]-asyoffset,Xup[i],Yup[j]);
			}else{
				pad[i][j] = new TPad(padName.Data(),padName.Data(),
									 Xlow[i],Ylow[j],Xup[i],Yup[j]-asyoffset);
			}
			
			
			if(i==0) pad[i][j]->SetLeftMargin(leftMargin);
			else pad[i][j]->SetLeftMargin(0);
			
			if(i==(columns-1)) pad[i][j]->SetRightMargin(edge);
			else pad[i][j]->SetRightMargin(0);
			
			if(j==0) pad[i][j]->SetTopMargin(edge);
			//else pad[i][j]->SetTopMargin(0.01);
			else pad[i][j]->SetTopMargin(0.02);
			
			//if(j==0) pad[i][j]->SetTopMargin(edge*3.5);
			//else pad[i][j]->SetTopMargin(0.0);
			
			if(j==(rows-1)) pad[i][j]->SetBottomMargin(bottomMargin);
			else pad[i][j]->SetBottomMargin(0.15);
			
			pad[i][j]->Draw();
			pad[i][j]->cd();
			pad[i][j]->SetNumber(columns*j+i+1);
			
		}
	}
}

void putCMSPrel(double x, double y, double size){
	TLatex *tex=0;
	tex = new TLatex(x,y,"CMS Preliminary");
	tex->SetTextSize(size);
	tex->SetLineWidth(2);
	tex->SetNDC();
	tex->Draw();
}
void drawText(const char *text, float xp, float yp, int size){
	TLatex *tex = new TLatex(xp,yp,text);
	tex->SetTextFont(63);
	tex->SetTextSize(size);
	tex->SetTextColor(kBlack);
	tex->SetLineWidth(1);
	//tex->SetTextFont(42);
	tex->SetNDC();
	tex->Draw();
}


void prepareNcollUnc(int nbins, float maxpt=300.){
	
	int fillsty = 1001;
	
	const int n = nbins;
	
	double xvalue[n];
	double yvalue[n];
	double xerror[n];
	double yerror1[n], yerror2[n], yerror3[n], yerror4[n], yerror5[n], yerror6[n];
	


  for(int i=0;i<nbins;i++){

    xvalue[i] = 300.1 + 1.2*(double)i, yvalue[i]=1.0, xerror[i]=0.0;  
    
      // TAA
    yerror1[i]=0.0409, yerror2[i]=0.0459, yerror3[i]=0.0578, yerror4[i]=0.0944, yerror5[i]=0.143, yerror6[i]=0.176;
 
    
    // add 6% error 
    yerror1[i]=TMath::Sqrt(yerror1[i]*yerror1[i]+0.06*0.06);
    yerror2[i]=TMath::Sqrt(yerror2[i]*yerror2[i]+0.06*0.06);
    yerror3[i]=TMath::Sqrt(yerror3[i]*yerror3[i]+0.06*0.06);
    yerror4[i]=TMath::Sqrt(yerror4[i]*yerror4[i]+0.06*0.06);
    yerror5[i]=TMath::Sqrt(yerror5[i]*yerror5[i]+0.06*0.06);
    yerror6[i]=TMath::Sqrt(yerror6[i]*yerror6[i]+0.06*0.06);
    cout<<"TAA + Lumi uncertainty = "<<yerror1[i]<<endl;   
  
  }
  
  
 // int ci = 29;
  int ci = 15;

  tTAAerr[0] = new TGraphErrors(n,xvalue,yvalue,xerror,yerror1);
  tTAAerr[0]->SetFillColor(ci);
  tTAAerr[0]->SetLineColor(ci);
  tTAAerr[0]->SetFillStyle(fillsty);

  tTAAerr[1] = new TGraphErrors(n,xvalue,yvalue,xerror,yerror2);
  tTAAerr[1]->SetFillColor(ci);
  tTAAerr[1]->SetFillStyle(fillsty);

  tTAAerr[2] = new TGraphErrors(n,xvalue,yvalue,xerror,yerror3);
  tTAAerr[2] ->SetFillColor(ci);
  tTAAerr[2] ->SetFillStyle(fillsty);

  tTAAerr[3]  = new TGraphErrors(n,xvalue,yvalue,xerror,yerror4);
  tTAAerr[3]->SetFillColor(ci);
  tTAAerr[3]->SetFillStyle(fillsty);

  tTAAerr[4] = new TGraphErrors(n,xvalue,yvalue,xerror,yerror5);
  tTAAerr[4]->SetFillColor(ci);
  tTAAerr[4]->SetFillStyle(fillsty);

  tTAAerr[5] = new TGraphErrors(n,xvalue,yvalue,xerror,yerror6);
  tTAAerr[5]->SetFillColor(ci);
  tTAAerr[5]->SetFillStyle(fillsty);
  


}
/*
void DrawNpartTAABand(){
	double xvalueNpart[6];
	double yerrorNpart[6];
	xvalueNpart[0] = 381.29; xvalueNpart[1] = 329.41; xvalueNpart[2] = 224.28;
  	xvalueNpart[3] = 108.12; xvalueNpart[4] = 42.04;  xvalueNpart[5] = 11.43;
  	yerrorNpart[0]=0.0409, yerrorNpart[1]=0.0459, yerrorNpart[2]=0.0578, yerrorNpart[3]=0.0944, yerrorNpart[4]=0.143, yerrorNpart[5]=0.176;
  	
  	int ci = 30;
  	
		for (int i=0;i<6;i++) {
		
			TBox *b = new TBox(xvalueNpart[i]-5,1.-yerrorNpart[i]/2,xvalueNpart[i]+5,1.+yerrorNpart[i]/2);
			b->SetFillColor(ci);
			b->SetFillStyle(3001);
			b->SetLineColor(ci);
			b->Draw();
		}
 
 
 }
 
*/
 

void dumpDatatoTxt(const char *centbin,TH1F *h, TH1F *hsys, TH1F *htotStat, const char *txtfile)
{
ofstream outf(txtfile,ios::out);
for(int ix=1;ix<=h->GetNbinsX();ix++){
double pt = h->GetBinCenter(ix);
double val = h->GetBinContent(ix);
double Uncorerr = h->GetBinError(ix);
double syserr = hsys->GetBinContent(ix)-1;
double totStaterr = htotStat->GetBinError(ix);

outf<<setprecision(0) << fixed <<pt<<"\t" << setprecision(3) << fixed <<val<<"\t" << setprecision(5) << fixed << totStaterr<<"\t"<< setprecision(4) << fixed << syserr*val << endl;
}
outf.close();
}


TGraphErrors *ShiftGraph(TGraphErrors* pGraph, Double_t pNumber){
    // shifts a graph by the absolute value in the argument
    TGraphErrors *pGraphtmp;
    for (Int_t i=0;i<pGraph->GetN();i++){
	Double_t x,y;
	Double_t yerr;
	pGraph->GetPoint(i,x,y);
	yerr = pGraph->GetErrorY(i);
	x = x+pNumber;
	pGraphtmp->SetPoint(i,x,y);
	pGraphtmp->SetPointError(i,yerr);
	}
	return pGraphtmp;
}


TGraphErrors* HistToTgraphShift(TH1F* hist,Double_t pNumber){

  TGraphErrors *pGraphtmp;
  int nbins = hist->GetNbinsX();

  const int nlines = nbins;

  float pt[nlines], xsec[nlines];
  float pterr[nlines], xsecerr[nlines];

  for(int i = 0; i<nbins; i++ ){
    pt[i] = hist->GetBinCenter(i+1);
    xsec[i] = hist->GetBinContent(i+1);
    xsecerr[i] = hist->GetBinError(i+1);
    pterr[i] = 0;
  }

  pGraphtmp = new TGraphErrors(nlines,pt,xsec,pterr,xsecerr);
  
  for (Int_t i=0;i<pGraphtmp->GetN();i++){
	Double_t x,y;
	Double_t yerr;
	pGraphtmp->GetPoint(i,x,y);
	yerr = pGraphtmp->GetErrorY(i);
	x = x+pNumber;
	pGraphtmp->SetPoint(i,x,y);
	pGraphtmp->SetPointError(i,yerr);
	}
  return pGraphtmp;
}


void DrawPanelLabel(int i){
 	TLatex *tex; 
	
 	 if(i==0) tex = new TLatex(0.07,0.9,"(f)");
	 if(i==1) tex = new TLatex(0.07,0.9,"(e)");
	 if(i==2) tex = new TLatex(0.19,0.9,"(d)");
	 if(i==3) tex = new TLatex(0.07,0.9,"(c)");
	 if(i==4) tex = new TLatex(0.07,0.9,"(b)");
	 if(i==5) tex = new TLatex(0.19,0.9,"(a)");
	 tex->SetTextFont(63);
	 tex->SetTextSize(18);
	 tex->SetTextColor(kBlack);
	 tex->SetLineWidth(1);
	 tex->SetNDC();

	 
 	tex->Draw();
	
 }
