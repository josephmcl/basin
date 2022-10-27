#include "hybrid_sbp_sat_2d.h"

std::vector<std::size_t> rowidx = {
  9, 10, 11, 122, 123, 124, 20, 21, 22, 133, 134, 135, 31, 32, 33, 144, 145, 146, 42, 43, 44, 155, 156, 157, 53, 54, 55, 166, 167, 168, 64, 65, 66, 177, 178, 179, 75, 76, 77, 188, 189, 190, 86, 87, 88, 199, 200, 201, 97, 98, 99, 210, 211, 212, 108, 109, 110, 221, 222, 223, 119, 120, 121, 232, 233, 234, 130, 131, 132, 243, 244, 245, 141, 142, 143, 254, 255, 256, 152, 153, 154, 265, 266, 267, 163, 164, 165, 276, 277, 278, 174, 175, 176, 287, 288, 289, 185, 186, 187, 298, 299, 300, 196, 197, 198, 309, 310, 311, 207, 208, 209, 320, 321, 322, 218, 219, 220, 331, 332, 333, 229, 230, 231, 342, 343, 344, 240, 241, 242, 353, 354, 355, 89, 100, 111, 364, 375, 386, 90, 101, 112, 365, 376, 387, 91, 102, 113, 366, 377, 388, 92, 103, 114, 367, 378, 389, 93, 104, 115, 368, 379, 390, 94, 105, 116, 369, 380, 391, 95, 106, 117, 370, 381, 392, 96, 107, 118, 371, 382, 393, 97, 108, 119, 372, 383, 394, 98, 109, 120, 373, 384, 395, 99, 110, 121, 374, 385, 396, 210, 221, 232, 485, 496, 507, 211, 222, 233, 486, 497, 508, 212, 223, 234, 487, 498, 509, 213, 224, 235, 488, 499, 510, 214, 225, 236, 489, 500, 511, 215, 226, 237, 490, 501, 512, 216, 227, 238, 491, 502, 513, 217, 228, 239, 492, 503, 514, 218, 229, 240, 493, 504, 515, 219, 230, 241, 494, 505, 516, 220, 231, 242, 495, 506, 517, 331, 342, 353, 606, 617, 628, 332, 343, 354, 607, 618, 629, 333, 344, 355, 608, 619, 630, 334, 345, 356, 609, 620, 631, 335, 346, 357, 610, 621, 632, 336, 347, 358, 611, 622, 633, 337, 348, 359, 612, 623, 634, 338, 349, 360, 613, 624, 635, 339, 350, 361, 614, 625, 636, 340, 351, 362, 615, 626, 637, 341, 352, 363, 616, 627, 638, 372, 373, 374, 485, 486, 487, 383, 384, 385, 496, 497, 498, 394, 395, 396, 507, 508, 509, 405, 406, 407, 518, 519, 520, 416, 417, 418, 529, 530, 531, 427, 428, 429, 540, 541, 542, 438, 439, 440, 551, 552, 553, 449, 450, 451, 562, 563, 564, 460, 461, 462, 573, 574, 575, 471, 472, 473, 584, 585, 586, 482, 483, 484, 595, 596, 597, 493, 494, 495, 606, 607, 608, 504, 505, 506, 617, 618, 619, 515, 516, 517, 628, 629, 630, 526, 527, 528, 639, 640, 641, 537, 538, 539, 650, 651, 652, 548, 549, 550, 661, 662, 663, 559, 560, 561, 672, 673, 674, 570, 571, 572, 683, 684, 685, 581, 582, 583, 694, 695, 696, 592, 593, 594, 705, 706, 707, 603, 604, 605, 716, 717, 718, 452, 463, 474, 727, 738, 749, 453, 464, 475, 728, 739, 750, 454, 465, 476, 729, 740, 751, 455, 466, 477, 730, 741, 752, 456, 467, 478, 731, 742, 753, 457, 468, 479, 732, 743, 754, 458, 469, 480, 733, 744, 755, 459, 470, 481, 734, 745, 756, 460, 471, 482, 735, 746, 757, 461, 472, 483, 736, 747, 758, 462, 473, 484, 737, 748, 759, 573, 584, 595, 848, 859, 870, 574, 585, 596, 849, 860, 871, 575, 586, 597, 850, 861, 872, 576, 587, 598, 851, 862, 873, 577, 588, 599, 852, 863, 874, 578, 589, 600, 853, 864, 875, 579, 590, 601, 854, 865, 876, 580, 591, 602, 855, 866, 877, 581, 592, 603, 856, 867, 878, 582, 593, 604, 857, 868, 879, 583, 594, 605, 858, 869, 880, 694, 705, 716, 969, 980, 991, 695, 706, 717, 970, 981, 992, 696, 707, 718, 971, 982, 993, 697, 708, 719, 972, 983, 994, 698, 709, 720, 973, 984, 995, 699, 710, 721, 974, 985, 996, 700, 711, 722, 975, 986, 997, 701, 712, 723, 976, 987, 998, 702, 713, 724, 977, 988, 999, 703, 714, 725, 978, 989, 1000, 704, 715, 726, 979, 990, 1001, 735, 736, 737, 848, 849, 850, 746, 747, 748, 859, 860, 861, 757, 758, 759, 870, 871, 872, 768, 769, 770, 881, 882, 883, 779, 780, 781, 892, 893, 894, 790, 791, 792, 903, 904, 905, 801, 802, 803, 914, 915, 916, 812, 813, 814, 925, 926, 927, 823, 824, 825, 936, 937, 938, 834, 835, 836, 947, 948, 949, 845, 846, 847, 958, 959, 960, 856, 857, 858, 969, 970, 971, 867, 868, 869, 980, 981, 982, 878, 879, 880, 991, 992, 993, 889, 890, 891, 1002, 1003, 1004, 900, 901, 902, 1013, 1014, 1015, 911, 912, 913, 1024, 1025, 1026, 922, 923, 924, 1035, 1036, 1037, 933, 934, 935, 1046, 1047, 1048, 944, 945, 946, 1057, 1058, 1059, 955, 956, 957, 1068, 1069, 1070, 966, 967, 968, 1079, 1080, 1081
};

std::vector<std::size_t> colidx = {
  1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 8, 8, 8, 8, 8, 8, 9, 9, 9, 9, 9, 9, 10, 10, 10, 10, 10, 10, 11, 11, 11, 11, 11, 11, 12, 12, 12, 12, 12, 12, 13, 13, 13, 13, 13, 13, 14, 14, 14, 14, 14, 14, 15, 15, 15, 15, 15, 15, 16, 16, 16, 16, 16, 16, 17, 17, 17, 17, 17, 17, 18, 18, 18, 18, 18, 18, 19, 19, 19, 19, 19, 19, 20, 20, 20, 20, 20, 20, 21, 21, 21, 21, 21, 21, 22, 22, 22, 22, 22, 22, 23, 23, 23, 23, 23, 23, 24, 24, 24, 24, 24, 24, 25, 25, 25, 25, 25, 25, 26, 26, 26, 26, 26, 26, 27, 27, 27, 27, 27, 27, 28, 28, 28, 28, 28, 28, 29, 29, 29, 29, 29, 29, 30, 30, 30, 30, 30, 30, 31, 31, 31, 31, 31, 31, 32, 32, 32, 32, 32, 32, 33, 33, 33, 33, 33, 33, 34, 34, 34, 34, 34, 34, 35, 35, 35, 35, 35, 35, 36, 36, 36, 36, 36, 36, 37, 37, 37, 37, 37, 37, 38, 38, 38, 38, 38, 38, 39, 39, 39, 39, 39, 39, 40, 40, 40, 40, 40, 40, 41, 41, 41, 41, 41, 41, 42, 42, 42, 42, 42, 42, 43, 43, 43, 43, 43, 43, 44, 44, 44, 44, 44, 44, 45, 45, 45, 45, 45, 45, 46, 46, 46, 46, 46, 46, 47, 47, 47, 47, 47, 47, 48, 48, 48, 48, 48, 48, 49, 49, 49, 49, 49, 49, 50, 50, 50, 50, 50, 50, 51, 51, 51, 51, 51, 51, 52, 52, 52, 52, 52, 52, 53, 53, 53, 53, 53, 53, 54, 54, 54, 54, 54, 54, 55, 55, 55, 55, 55, 55, 56, 56, 56, 56, 56, 56, 57, 57, 57, 57, 57, 57, 58, 58, 58, 58, 58, 58, 59, 59, 59, 59, 59, 59, 60, 60, 60, 60, 60, 60, 61, 61, 61, 61, 61, 61, 62, 62, 62, 62, 62, 62, 63, 63, 63, 63, 63, 63, 64, 64, 64, 64, 64, 64, 65, 65, 65, 65, 65, 65, 66, 66, 66, 66, 66, 66, 67, 67, 67, 67, 67, 67, 68, 68, 68, 68, 68, 68, 69, 69, 69, 69, 69, 69, 70, 70, 70, 70, 70, 70, 71, 71, 71, 71, 71, 71, 72, 72, 72, 72, 72, 72, 73, 73, 73, 73, 73, 73, 74, 74, 74, 74, 74, 74, 75, 75, 75, 75, 75, 75, 76, 76, 76, 76, 76, 76, 77, 77, 77, 77, 77, 77, 78, 78, 78, 78, 78, 78, 79, 79, 79, 79, 79, 79, 80, 80, 80, 80, 80, 80, 81, 81, 81, 81, 81, 81, 82, 82, 82, 82, 82, 82, 83, 83, 83, 83, 83, 83, 84, 84, 84, 84, 84, 84, 85, 85, 85, 85, 85, 85, 86, 86, 86, 86, 86, 86, 87, 87, 87, 87, 87, 87, 88, 88, 88, 88, 88, 88, 89, 89, 89, 89, 89, 89, 90, 90, 90, 90, 90, 90, 91, 91, 91, 91, 91, 91, 92, 92, 92, 92, 92, 92, 93, 93, 93, 93, 93, 93, 94, 94, 94, 94, 94, 94, 95, 95, 95, 95, 95, 95, 96, 96, 96, 96, 96, 96, 97, 97, 97, 97, 97, 97, 98, 98, 98, 98, 98, 98, 99, 99, 99, 99, 99, 99, 100, 100, 100, 100, 100, 100, 101, 101, 101, 101, 101, 101, 102, 102, 102, 102, 102, 102, 103, 103, 103, 103, 103, 103, 104, 104, 104, 104, 104, 104, 105, 105, 105, 105, 105, 105, 106, 106, 106, 106, 106, 106, 107, 107, 107, 107, 107, 107, 108, 108, 108, 108, 108, 108, 109, 109, 109, 109, 109, 109, 110, 110, 110, 110, 110, 110, 111, 111, 111, 111, 111, 111, 112, 112, 112, 112, 112, 112, 113, 113, 113, 113, 113, 113, 114, 114, 114, 114, 114, 114, 115, 115, 115, 115, 115, 115, 116, 116, 116, 116, 116, 116, 117, 117, 117, 117, 117, 117, 118, 118, 118, 118, 118, 118, 119, 119, 119, 119, 119, 119, 120, 120, 120, 120, 120, 120, 121, 121, 121, 121, 121, 121, 122, 122, 122, 122, 122, 122, 123, 123, 123, 123, 123, 123, 124, 124, 124, 124, 124, 124, 125, 125, 125, 125, 125, 125, 126, 126, 126, 126, 126, 126, 127, 127, 127, 127, 127, 127, 128, 128, 128, 128, 128, 128, 129, 129, 129, 129, 129, 129, 130, 130, 130, 130, 130, 130, 131, 131, 131, 131, 131, 131, 132, 132, 132, 132, 132, 132
  };

std::vector<long double> ftvals = {
  0.25, -1.0, 0.7333333333333333, 0.7333333333333333, -1.0, 0.25, 0.5, 
  -2.0, 1.4666666666666666, 1.4666666666666666, -2.0, 0.5, 0.5, -2.0, 1.4666666666666666, 1.4666666666666666, -2.0, 0.5, 0.5, -2.0, 1.4666666666666666, 1.4666666666666666, -2.0, 0.5, 0.5, -2.0, 1.4666666666666666, 1.4666666666666666, -2.0, 0.5, 0.5, -2.0, 1.4666666666666666, 1.4666666666666666, -2.0, 0.5, 0.5, -2.0, 1.4666666666666666, 1.4666666666666666, -2.0, 0.5, 0.5, -2.0, 1.4666666666666666, 1.4666666666666666, -2.0, 0.5, 0.5, -2.0, 1.4666666666666666, 1.4666666666666666, -2.0, 0.5, 0.5, -2.0, 1.4666666666666666, 1.4666666666666666, -2.0, 0.5, 0.25, -1.0, 0.7333333333333333, 0.7333333333333333, -1.0, 0.25, 0.25, -1.0, 0.7333333333333333, 0.7333333333333333, -1.0, 0.25, 0.5, -2.0, 1.4666666666666666, 1.4666666666666666, -2.0, 0.5, 0.5, -2.0, 1.4666666666666666, 1.4666666666666666, -2.0, 0.5, 0.5, -2.0, 1.4666666666666666, 1.4666666666666666, -2.0, 0.5, 0.5, -2.0, 1.4666666666666666, 1.4666666666666666, -2.0, 0.5, 0.5, -2.0, 1.4666666666666666, 1.4666666666666666, -2.0, 0.5, 0.5, -2.0, 1.4666666666666666, 1.4666666666666666, -2.0, 0.5, 0.5, -2.0, 1.4666666666666666, 1.4666666666666666, -2.0, 0.5, 0.5, -2.0, 1.4666666666666666, 1.4666666666666666, -2.0, 0.5, 0.5, -2.0, 1.4666666666666666, 1.4666666666666666, -2.0, 0.5, 0.25, -1.0, 0.7333333333333333, 0.7333333333333333, -1.0, 0.25, 0.25, -1.0, 0.7333333333333333, 0.7333333333333333, -1.0, 0.25, 0.5, -2.0, 1.4666666666666666, 1.4666666666666666, -2.0, 0.5, 0.5, -2.0, 1.4666666666666666, 1.4666666666666666, -2.0, 0.5, 0.5, -2.0, 1.4666666666666666, 1.4666666666666666, -2.0, 0.5, 0.5, -2.0, 1.4666666666666666, 1.4666666666666666, -2.0, 0.5, 0.5, -2.0, 1.4666666666666666, 1.4666666666666666, -2.0, 0.5, 0.5, -2.0, 1.4666666666666666, 1.4666666666666666, -2.0, 0.5, 0.5, -2.0, 1.4666666666666666, 1.4666666666666666, -2.0, 0.5, 0.5, -2.0, 1.4666666666666666, 1.4666666666666666, -2.0, 0.5, 0.5, -2.0, 1.4666666666666666, 1.4666666666666666, -2.0, 0.5, 0.25, -1.0, 0.7333333333333333, 0.7333333333333333, -1.0, 0.25, 0.25, -1.0, 0.7333333333333333, 0.7333333333333333, -1.0, 0.25, 0.5, -2.0, 1.4666666666666666, 1.4666666666666666, -2.0, 0.5, 0.5, -2.0, 1.4666666666666666, 1.4666666666666666, -2.0, 0.5, 0.5, -2.0, 1.4666666666666666, 1.4666666666666666, -2.0, 0.5, 0.5, -2.0, 1.4666666666666666, 1.4666666666666666, -2.0, 0.5, 0.5, -2.0, 1.4666666666666666, 1.4666666666666666, -2.0, 0.5, 0.5, -2.0, 1.4666666666666666, 1.4666666666666666, -2.0, 0.5, 0.5, -2.0, 1.4666666666666666, 1.4666666666666666, -2.0, 0.5, 0.5, -2.0, 1.4666666666666666, 1.4666666666666666, -2.0, 0.5, 0.5, -2.0, 1.4666666666666666, 1.4666666666666666, -2.0, 0.5, 0.25, -1.0, 0.7333333333333333, 0.7333333333333333, -1.0, 0.25, 0.25, -1.0, 0.7333333333333333, 0.7333333333333333, -1.0, 0.25, 0.5, -2.0, 1.4666666666666666, 1.4666666666666666, -2.0, 0.5, 0.5, -2.0, 1.4666666666666666, 1.4666666666666666, -2.0, 0.5, 0.5, -2.0, 1.4666666666666666, 1.4666666666666666, -2.0, 0.5, 0.5, -2.0, 1.4666666666666666, 1.4666666666666666, -2.0, 0.5, 0.5, -2.0, 1.4666666666666666, 1.4666666666666666, -2.0, 0.5, 0.5, -2.0, 1.4666666666666666, 1.4666666666666666, -2.0, 0.5, 0.5, -2.0, 1.4666666666666666, 1.4666666666666666, -2.0, 0.5, 0.5, -2.0, 1.4666666666666666, 1.4666666666666666, -2.0, 0.5, 0.5, -2.0, 1.4666666666666666, 1.4666666666666666, -2.0, 0.5, 0.25, -1.0, 0.7333333333333333, 0.7333333333333333, -1.0, 0.25, 0.25, -1.0, 0.7333333333333333, 0.7333333333333333, -1.0, 0.25, 0.5, -2.0, 1.4666666666666666, 1.4666666666666666, -2.0, 0.5, 0.5, -2.0, 1.4666666666666666, 1.4666666666666666, -2.0, 0.5, 0.5, -2.0, 1.4666666666666666, 1.4666666666666666, -2.0, 0.5, 0.5, -2.0, 1.4666666666666666, 1.4666666666666666, -2.0, 0.5, 0.5, -2.0, 1.4666666666666666, 1.4666666666666666, -2.0, 0.5, 0.5, -2.0, 1.4666666666666666, 1.4666666666666666, -2.0, 0.5, 0.5, -2.0, 1.4666666666666666, 1.4666666666666666, -2.0, 0.5, 0.5, -2.0, 1.4666666666666666, 1.4666666666666666, -2.0, 0.5, 0.5, -2.0, 1.4666666666666666, 1.4666666666666666, -2.0, 0.5, 0.25, -1.0, 0.7333333333333333, 0.7333333333333333, -1.0, 0.25, 0.25, -1.0, 0.7333333333333333, 0.7333333333333333, -1.0, 0.25, 0.5, -2.0, 1.4666666666666666, 1.4666666666666666, -2.0, 0.5, 0.5, -2.0, 1.4666666666666666, 1.4666666666666666, -2.0, 0.5, 0.5, -2.0, 1.4666666666666666, 1.4666666666666666, -2.0, 0.5, 0.5, -2.0, 1.4666666666666666, 1.4666666666666666, -2.0, 0.5, 0.5, -2.0, 1.4666666666666666, 1.4666666666666666, -2.0, 0.5, 0.5, -2.0, 1.4666666666666666, 1.4666666666666666, -2.0, 0.5, 0.5, -2.0, 1.4666666666666666, 1.4666666666666666, -2.0, 0.5, 0.5, -2.0, 1.4666666666666666, 1.4666666666666666, -2.0, 0.5, 0.5, -2.0, 1.4666666666666666, 1.4666666666666666, -2.0, 0.5, 0.25, -1.0, 0.7333333333333333, 0.7333333333333333, -1.0, 0.25, 0.25, -1.0, 0.7333333333333333, 0.7333333333333333, -1.0, 0.25, 0.5, -2.0, 1.4666666666666666, 1.4666666666666666, -2.0, 0.5, 0.5, -2.0, 1.4666666666666666, 1.4666666666666666, -2.0, 0.5, 0.5, -2.0, 1.4666666666666666, 1.4666666666666666, -2.0, 0.5, 0.5, -2.0, 1.4666666666666666, 1.4666666666666666, -2.0, 0.5, 0.5, -2.0, 1.4666666666666666, 1.4666666666666666, -2.0, 0.5, 0.5, -2.0, 1.4666666666666666, 1.4666666666666666, -2.0, 0.5, 0.5, -2.0, 1.4666666666666666, 1.4666666666666666, -2.0, 0.5, 0.5, -2.0, 1.4666666666666666, 1.4666666666666666, -2.0, 0.5, 0.5, -2.0, 1.4666666666666666, 1.4666666666666666, -2.0, 0.5, 0.25, -1.0, 0.7333333333333333, 0.7333333333333333, -1.0, 0.25, 0.25, -1.0, 0.7333333333333333, 0.7333333333333333, -1.0, 0.25, 0.5, -2.0, 1.4666666666666666, 1.4666666666666666, -2.0, 0.5, 0.5, -2.0, 1.4666666666666666, 1.4666666666666666, -2.0, 0.5, 0.5, -2.0, 1.4666666666666666, 1.4666666666666666, -2.0, 0.5, 0.5, -2.0, 1.4666666666666666, 1.4666666666666666, -2.0, 0.5, 0.5, -2.0, 1.4666666666666666, 1.4666666666666666, -2.0, 0.5, 0.5, -2.0, 1.4666666666666666, 1.4666666666666666, -2.0, 0.5, 0.5, -2.0, 1.4666666666666666, 1.4666666666666666, -2.0, 0.5, 0.5, -2.0, 1.4666666666666666, 1.4666666666666666, -2.0, 0.5, 0.5, -2.0, 1.4666666666666666, 1.4666666666666666, -2.0, 0.5, 0.25, -1.0, 0.7333333333333333, 0.7333333333333333, -1.0, 0.25, 0.25, -1.0, 0.7333333333333333, 0.7333333333333333, -1.0, 0.25, 0.5, -2.0, 1.4666666666666666, 1.4666666666666666, -2.0, 0.5, 0.5, -2.0, 1.4666666666666666, 1.4666666666666666, -2.0, 0.5, 0.5, -2.0, 1.4666666666666666, 1.4666666666666666, -2.0, 0.5, 0.5, -2.0, 1.4666666666666666, 1.4666666666666666, -2.0, 0.5, 0.5, -2.0, 1.4666666666666666, 1.4666666666666666, -2.0, 0.5, 0.5, -2.0, 1.4666666666666666, 1.4666666666666666, -2.0, 0.5, 0.5, -2.0, 1.4666666666666666, 1.4666666666666666, -2.0, 0.5, 0.5, -2.0, 1.4666666666666666, 1.4666666666666666, -2.0, 0.5, 0.5, -2.0, 1.4666666666666666, 1.4666666666666666, -2.0, 0.5, 0.25, -1.0, 0.7333333333333333, 0.7333333333333333, -1.0, 0.25, 0.25, -1.0, 0.7333333333333333, 0.7333333333333333, -1.0, 0.25, 0.5, -2.0, 1.4666666666666666, 1.4666666666666666, -2.0, 0.5, 0.5, -2.0, 1.4666666666666666, 1.4666666666666666, -2.0, 0.5, 0.5, -2.0, 1.4666666666666666, 1.4666666666666666, -2.0, 0.5, 0.5, -2.0, 1.4666666666666666, 1.4666666666666666, -2.0, 0.5, 0.5, -2.0, 1.4666666666666666, 1.4666666666666666, -2.0, 0.5, 0.5, -2.0, 1.4666666666666666, 1.4666666666666666, -2.0, 0.5, 0.5, -2.0, 1.4666666666666666, 1.4666666666666666, -2.0, 0.5, 0.5, -2.0, 1.4666666666666666, 1.4666666666666666, -2.0, 0.5, 0.5, -2.0, 1.4666666666666666, 1.4666666666666666, -2.0, 0.5, 0.25, -1.0, 0.7333333333333333, 0.7333333333333333, -1.0, 0.25, 0.25, -1.0, 0.7333333333333333, 0.7333333333333333, -1.0, 0.25, 0.5, -2.0, 1.4666666666666666, 1.4666666666666666, -2.0, 0.5, 0.5, -2.0, 1.4666666666666666, 1.4666666666666666, -2.0, 0.5, 0.5, -2.0, 1.4666666666666666, 1.4666666666666666, -2.0, 0.5, 0.5, -2.0, 1.4666666666666666, 1.4666666666666666, -2.0, 0.5, 0.5, -2.0, 1.4666666666666666, 1.4666666666666666, -2.0, 0.5, 0.5, -2.0, 1.4666666666666666, 1.4666666666666666, -2.0, 0.5, 0.5, -2.0, 1.4666666666666666, 1.4666666666666666, -2.0, 0.5, 0.5, -2.0, 1.4666666666666666, 1.4666666666666666, -2.0, 0.5, 0.5, -2.0, 1.4666666666666666, 1.4666666666666666, -2.0, 0.5, 0.25, -1.0, 0.7333333333333333, 0.7333333333333333, -1.0, 0.25};

/* Given three petsc matrices, return three portions of the F matrix for 
   a fixed problem size. TODO: generalize this solution for and B and any
   N, i.e., number of blocks and global problem size. */
void lazy_f(
  sbp_sat::petsc_matrix &f1,
  sbp_sat::petsc_matrix &f2,
  sbp_sat::petsc_matrix &f3) { //1089


    linalg::make_local_sparse_matrix<sbp_sat::x2::fw>(f1, 363, 132, 4);
    linalg::make_local_sparse_matrix<sbp_sat::x2::fw>(f2, 363, 132, 4);
    linalg::make_local_sparse_matrix<sbp_sat::x2::fw>(f3, 363, 132, 4);
    for (std::size_t i = 0; i < rowidx.size(); ++i) {

      if (rowidx[i] < 363) {
        linalg::set_matrix_value<sbp_sat::x2::fw>(
          f1, rowidx[i] - 1, colidx[i] - 1, ftvals[i]);
      }
      else if (rowidx[i] < 726) {
        linalg::set_matrix_value<sbp_sat::x2::fw>(
          f1, rowidx[i] - 364, colidx[i] - 1, ftvals[i]);
      }
      else {
        linalg::set_matrix_value<sbp_sat::x2::fw>(
          f3, rowidx[i] - 727, colidx[i] - 1, ftvals[i]);
      }
    }
    linalg::finalize<sbp_sat::x2::fw>(f1);
    linalg::finalize<sbp_sat::x2::fw>(f2);
    linalg::finalize<sbp_sat::x2::fw>(f3);
}

void lazy_ft(
  std::vector<sbp_sat::petsc_matrix> &f,
  std::vector<sbp_sat::petsc_matrix> &ft) { //1089

    for (auto i = std::size_t(0); i != f.size(); ++i) {
      MatTranspose(f[i], MAT_INITIAL_MATRIX, &ft[ft.size() - i]);
      linalg::finalize<sbp_sat::x2::fw>(ft[ft.size() - i]);
    }
}

/* [ M  F ][ u ] = [ g   ]
     [ FT D ][ λ ]   [ g_λ ]
  
    [ B   T ][ u ] = [ bu ]
    [ TT  D ][ λ ]   [ bλ ]
  */

/* Baseline Poisson solution. Solves A u = b given a 2D volume point 
   dimension and boundary conditions. */
void sbp_sat::x2::
petsc_poisson(sbp_sat::real_v             &result,
              sbp_sat::domain_v     const &domain, 
              sbp_sat::block_t      const &block, 
              sbp_sat::boundary_tx2 const &boundaries) {



   
  using namespace sbp_sat;

  auto A = matrix<fw> {};

  auto w = neumann;
  auto e = neumann;
  auto n = dirichlet;
  auto s = dirichlet;

  x2::write_m(A, block, block, {w, e, n, s});

  finalize<fw>(A);

  destroy<fw>(A);

}

void sbp_sat::x2::
petsc_hybridized_poisson(sbp_sat::real_v             &result,
                         sbp_sat::domain_v     const &domain, 
                         sbp_sat::block_v      const &blocks, 
                         sbp_sat::boundary_vx2 const &boundaries) {


  /* [ M  F ][ u ] = [ g   ]
     [ FT D ][ λ ]   [ g_λ ]
  
    [ B   T ][ u ] = [ bu ]
    [ TT  D ][ λ ]   [ bλ ]
  */
   
  using namespace sbp_sat;

  nat_t local_problems = blocks.size() * blocks.size();

  // only for square hybrid systems 
  const std::size_t n_blocks = blocks.size() * blocks.size();

  /*  Create petsc objects for the block diagonal matrix, M. 

          [ m_0           ] 
      M = [     ...       ] in ( M  F )  
          [         m_n-1 ]    ( Ft D )   

      From the segmented volume point blocks

      | b0  | ... | ...  |
      |-----|-----|------|
      | ... | ... | ...  |
      |-----|-----|------|
      | ... | ... | bn-1 |  */

  std::vector<std::vector<std::size_t>> interfaces = {
    {0, 1, 0, 7, 0, 0, 0, 0, 0}, 
    {0, 0, 2, 0, 8, 0, 0, 0, 0}, 
    {0, 0, 0, 0, 0, 9, 0, 0, 0}, 
    {0, 0, 0, 0, 3, 0, 10, 0, 0}, 
    {0, 0, 0, 0, 0, 4, 0, 11, 0}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 12}, 
    {0, 0, 0, 0, 0, 0, 0, 5, 0}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 6}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 0}};

  std::size_t n_interfaces = 0;
  for (auto &row : interfaces) {
    for (auto &element : row) {
      if (element != 0) {
        n_interfaces += 1;
      }
    }
  }

  constexpr std::size_t w = 1; constexpr std::size_t e = 2; 
  constexpr std::size_t s = 3; constexpr std::size_t n = 4; 

  std::vector<std::vector<std::size_t>> F_symbols(n_blocks, // rows 
    std::vector<std::size_t>(n_interfaces, 0)); // columns
  std::vector<std::vector<std::size_t>> FT_symbols(n_interfaces, // rows
    std::vector<std::size_t>(n_blocks, 0)); // columns

  // set F_symbols, FT_symbols
  for (std::size_t row = 0; row < n_blocks; ++row) {
    for (std::size_t col = 0; col < n_blocks; ++col) {
      std::size_t interface = interfaces[row][col];
      if (interface != 0 and row == col - 1) {  
        F_symbols[row][interface - 1] = n;
        F_symbols[col][interface - 1] = s;
        FT_symbols[interface - 1][row] = n;
        FT_symbols[interface - 1][col] = s;
      }
      else if (interface != 0 and row == col - 3) {  
        F_symbols[row][interface - 1] = w;
        F_symbols[col][interface - 1] = e;
        FT_symbols[interface - 1][row] = w;
        FT_symbols[interface - 1][col] = e;
      }
    }
  }

  auto [rows, cols] = blocks[0];
  auto spacing = rows.to - rows.from;
  auto span = spacing;
  auto n_points_x = rows.size();

  std::cout << "Square hybrid specs:" << std::endl 
    << " | local problem size: " << n_points_x << " x " << n_points_x 
    << std::endl << " | span: " << span 
    << std::endl << " | total blocks: " << n_blocks 
    << std::endl << " | total interfaces: " << n_interfaces 
    << std::endl;


  std::cout << n_points_x << ", " << span << std::endl;

  // n_points = number of columns in 
  //            matrix component. 

  // components class contains most of the sbp-sat component matrices
  // needed to set up an sbp-sat problem.
  auto sbp = components{11, span};

  sbp.τ = 42.; // hard code these coeffs for now. 
  sbp.β = 1.;

  // std::size_t g_size = n_blocks * n_interfaces;

  // F components are stored as vectors because we compute
  // M^-1 F by solving every row each component. 
  auto f = std::vector<std::vector<petsc_vector>>(
    4, std::vector<petsc_vector>(n_blocks));  
 
  make_f_subs(sbp, f);

  auto m = std::vector<petsc_matrix>(n_blocks);

  
  std::size_t bw, be, bn, bs;
  std::cout << local_problems << " local problems." << std::endl;
  for (std::size_t row = 0; row != blocks.size(); ++row) {
    for (std::size_t col = 0; col != blocks.size(); ++col) {

      auto matrix_index = blocks.size() * row + col;

      auto [row_x1, col_x1] = blocks[row];
      auto [row_x2, col_x2] = blocks[col]; 

      (void) col_x1; (void) row_x2; // unused

      bw = neumann;
      be = neumann;
      bn = col == 0                 ? dirichlet : neumann;
      bs = col == blocks.size() - 1 ? dirichlet : neumann;

      x2::write_m(m[matrix_index], blocks[row], blocks[col], {bw, be, bn, bs});

      auto m_index = (row * blocks.size()) + col + 1;
      auto ordinal = 
        m_index == 1 ? "st" : 
        m_index == 2 ? "nd" : 
        m_index == 3 ? "rd" : "th";

      std::cout << "Assembled " << m_index << ordinal << " M block (" 
        << row << ", " << col << ")" << std::endl;

    }
  }

  // MatView(m[n_blocks - 1], PETSC_VIEWER_STDOUT_SELF);
  
  
  std::cout << m.size() << " = m size " << std::endl;

  petsc_matrix II;
  MatCreateConstantDiagonal(PETSC_COMM_SELF, 
    11, 11, PETSC_DECIDE, PETSC_DECIDE, 1., &II);
  finalize<fw>(II);
  
  MatView(II, PETSC_VIEWER_STDOUT_SELF);

  // auto solvers = std::vector<KSP>(m.size());
  KSP solvers[9];
  for (std::size_t i = 0; i != m.size(); ++i) {
    solvers[i] = KSP();
    KSPCreate(PETSC_COMM_SELF, &solvers[i]); 
    KSPSetOperators(solvers[i], m[i], m[i]);
    // KSPSetUp(solvers[i]);
    int msz, nsz;
    MatGetSize(m[i], &msz, &nsz); std::cout 
    << "m size " << msz << ",  " << nsz << std::endl;
  }

  // Compute M^-1 F
  auto mf = std::vector<std::vector<petsc_vector>>(
    4 * n_blocks, std::vector<petsc_vector>(sbp.n)); 
  std::cout << "allocate mf" << std::endl; 
  for (std::size_t i = 0; i != 4 * n_blocks; ++i) {
    for (std::size_t j = 0; j != sbp.n; ++j) {
      VecCreateSeq(PETSC_COMM_SELF, sbp.n * sbp.n, &mf[i][j]);
    }
  }

  std::cout << "m size: " << m.size() << std::endl;
  std::cout << "f size: " << f.size() << ", " << f[0].size()<< std::endl;
  std::cout << "mf size: " << mf.size() << ", " << mf[0].size() << std::endl;

  msolvef(mf, &solvers[0], n_blocks, f);
  
  std::cout << "solved mf" << std::endl;

  std::cout << "end." << std::endl;
  for (std::size_t i = 0; i != 9; ++i) KSPDestroy(&solvers[i]);
  for (auto &e : m)  destroy<fw>(e);

  // STOP STOP STOP STOP STOP STOP STOP STOP STOP STOP STOP STOP 
  // STOP STOP STOP STOP STOP STOP STOP STOP STOP STOP STOP STOP 
  // STOP STOP STOP STOP STOP STOP STOP STOP STOP STOP STOP STOP 
  // STOP STOP STOP STOP STOP STOP STOP STOP STOP STOP STOP STOP 

  #ifdef SOMETHINGUNDEFINED

  auto M = std::vector<petsc_matrix>(local_problems); 

  

  auto F  = std::vector<petsc_matrix>(blocks.size()); 
  auto Ft = std::vector<petsc_matrix>(blocks.size()); 
  lazy_f(F[0], F[1], F[2]);
  std::cout << "Assembled F blocks (lazy.)" << std::endl;
  lazy_ft(F, Ft);
  std::cout << "Assembled Ft blocks (lazy.)" << std::endl;


  /* Set and finalize the g_bar vectors. */
  auto g_bar = std::vector<petsc_vector>(local_problems); 

  std::size_t g_size = blocks.size() * boundaries.size();

  for (std::size_t row = 0; row != blocks.size(); ++row) {
    for (std::size_t col = 0; col != blocks.size(); ++col) {
      // solve_d
    }
  }

  // PetscViewerPushFormat(PETSC_VIEWER_STD OUT_SELF, PETSC_VIEWER_ASCII_MATLAB);
  // MatView(a, PETSC_VIEWER_STDOUT_SELF);

  /* Free all petsc objects. */
  for (auto &e : M)  MatDestroy(&e);
  for (auto &e : F)  MatDestroy(&e);
  for (auto &e : Ft) MatDestroy(&e);
  for (auto &e : L) MatDestroy(&e);

  #endif
}

void sbp_sat::x2::
write_m(
  petsc_matrix       &m,
  block_t      const &block_x1, 
  block_t      const &block_x2, 
  std::array<std::size_t, 4> const bc) {

  auto [rows, c_] = block_x1;
  auto [r_, cols] = block_x2;
  (void) c_; (void) r_; // unused;

  auto x1_spacing = rows.to - rows.from;
  auto x2_spacing = cols.to - cols.from;

  auto x1_sz = rows.size() - 1;
  auto x2_sz = cols.size() - 1;

  auto x1_spacing_square = x1_spacing * x1_spacing;
  auto x2_spacing_square = x2_spacing * x2_spacing;


  auto hx1 = numerical::operators::H(rows.size(), 2, 0, x1_spacing);
  auto hx2 = numerical::operators::H(cols.size(), 2, 0, x2_spacing);

  // hard code for now 
  real_t τ = 1.;
  real_t β = 1.;

  // hard code the second order bs matrix for now. 
  vector_t bsx1 = {
    (3./2.) / x1_spacing * (x1_sz), 
        -2. / x1_spacing * (x1_sz),   // NOTE: maybe flip sz's? 
    (1./2.) / x1_spacing * (x1_sz)};

  vector_t bsx2 = {
    (3./2.) / x2_spacing * (x2_sz), 
        -2. / x2_spacing * (x2_sz), 
    (1./2.) / x2_spacing * (x2_sz)};

  matrix<fw> d2x1, d2x2, h1x1, h1x2, Ax1, Ax2;
  MatCreateSeqAIJ(PETSC_COMM_SELF, rows.size(), rows.size(), 4, 
    nullptr, &d2x1);
  MatCreateSeqAIJ(PETSC_COMM_SELF, cols.size(), cols.size(), 4, 
    nullptr, &d2x2);
  MatCreateSeqAIJ(PETSC_COMM_SELF, rows.size(), rows.size(), 1, 
    nullptr, &h1x1);
  MatCreateSeqAIJ(PETSC_COMM_SELF, cols.size(), cols.size(), 1, 
    nullptr, &h1x2);

  write_d2_h1(d2x1, hx1, rows, x1_spacing_square / (x1_sz * x1_sz), -1.);
  write_d2_h1(d2x2, hx2, cols, x2_spacing_square / (x2_sz * x2_sz), -1.);
  write_h1(h1x1, hx1);
  write_h1(h1x2, hx2);
  finalize<fw>(d2x1);
  finalize<fw>(d2x2);
  finalize<fw>(h1x1);
  finalize<fw>(h1x2);

  MatSeqAIJKron(h1x2, d2x1, MAT_INITIAL_MATRIX, &Ax1);
  MatSeqAIJKron(d2x2, h1x1, MAT_INITIAL_MATRIX, &Ax2);

  finalize<fw>(Ax2);
  finalize<fw>(Ax1);
  
  auto sz = rows.size() * cols.size();
  MatCreate(PETSC_COMM_SELF, &m);
  MatSetSizes(m, PETSC_DECIDE, PETSC_DECIDE, sz, sz);
  MatSetType(m, MATCOMPOSITE);

  MatCompositeAddMat(m, Ax1);
  MatCompositeAddMat(m, Ax2);

  // Append boundary condition coefficients as new matrices to the 
  // composite matrix, m. 
  add_boundary<x, left> (m, rows, cols, bsx1, hx2, β, τ, bc[0]);
  add_boundary<x, right>(m, rows, cols, bsx1, hx2, β, τ, bc[1]);
  add_boundary<y, left> (m, rows, cols, bsx1, hx2, β, τ, bc[2]);
  add_boundary<y, right>(m, rows, cols, bsx1, hx2, β, τ, bc[3]);
  
  MatCompositeSetMatStructure(m, DIFFERENT_NONZERO_PATTERN);
  MatCompositeMerge(m);
  finalize<fw>(m);

  destroy<fw>(d2x1);
  destroy<fw>(d2x2);
  destroy<fw>(h1x1);
  destroy<fw>(h1x2);
  destroy<fw>(Ax1);
  destroy<fw>(Ax2);
}

void sbp_sat::x2::write_d2_h1(
  petsc_matrix       &M, 
  real_v       const &h, 
  range_t      const &local,
  real_t       const  spacing_square, 
  real_t       const  coeff) {
    
  /* Initialize the first local skew row. */
  MatSetValue(M, 0, 0,  1. / spacing_square * h[0] * coeff, ADD_VALUES);
  MatSetValue(M, 0, 1, -2. / spacing_square * h[0] * coeff, ADD_VALUES);
  MatSetValue(M, 0, 2,  1. / spacing_square * h[0] * coeff, ADD_VALUES); 

  /* Initialize the final local skew row. */
  auto n = local.size() - 1;
  MatSetValue(M, n, n - 2,  1. / spacing_square * h[n] * coeff, ADD_VALUES);
  MatSetValue(M, n, n - 1, -2. / spacing_square * h[n] * coeff, ADD_VALUES);
  MatSetValue(M, n, n,      1. / spacing_square * h[n] * coeff, ADD_VALUES); 

  /* Initialize the interior local diagonal rows. */
  for (auto it = local.begin() + 1; it != local.end() - 1; ++it) {  
    auto i = it.index;
    MatSetValue(M, i, i - 1,  1. / spacing_square * h[i] * coeff, ADD_VALUES);
    MatSetValue(M, i, i,     -2. / spacing_square * h[i] * coeff, ADD_VALUES);
    MatSetValue(M, i, i + 1,  1. / spacing_square * h[i] * coeff, ADD_VALUES);  
  }
}

void sbp_sat::x2::write_h1(
  petsc_matrix                   &M, 
  std::vector<long double> const &h) {

  for (std::size_t i = 0; i < h.size(); ++i) {
    MatSetValue(M, i, i,  h[i], ADD_VALUES);
  }
}

void sbp_sat::x2::write_Ls(
  std::vector<sbp_sat::petsc_matrix> &L,
  sbp_sat::boundary_vx2 const &boundaries){
  for (std::size_t i = 0; i != boundaries.size(); ++i) {
    for (std::size_t j = 0; j != boundaries[i].size(); ++j) {
      auto [g_t, _] = boundaries[i][j];
      auto ell = &L[i + j * boundaries.size()];
      std::cout << i + j * boundaries.size() << std::endl;
      MatCreateSeqAIJ(PETSC_COMM_SELF, g_t.size(), 
        g_t.size() * g_t.size(), 4, nullptr, ell);

      if (j == 0) {
        for (std::size_t k = 0; k != g_t.size(); ++k)
          MatSetValue(*ell, k, k,  1., ADD_VALUES);
      }
      else if (j == 1) {
        for (std::size_t k = 0; k != g_t.size(); ++k)
          MatSetValue(*ell, k, g_t.size() * (g_t.size() - 1) + k,  1., ADD_VALUES);
      }
      else if (j == 2) {
        for (std::size_t k = 0; k != g_t.size(); ++k)
          MatSetValue(*ell, k, k * g_t.size(),  1., ADD_VALUES);
      }
      else if (j == 3) {
        for (std::size_t k = 0; k != g_t.size(); ++k)
          MatSetValue(*ell, k, (k + 1) * g_t.size() - 1,  1., ADD_VALUES);
      }
      linalg::finalize<sbp_sat::x2::fw>(*ell);
    }   
  }
}

void write_Lts(
  std::vector<sbp_sat::petsc_matrix> &L,
  std::vector<sbp_sat::petsc_matrix> &Lt) { 

    for (auto i = std::size_t(0); i != L.size(); ++i) {
      MatTranspose(L[i], MAT_INITIAL_MATRIX, &Lt[Lt.size() - i]);
      linalg::finalize<sbp_sat::x2::fw>(Lt[Lt.size() - i]);
    }
}

// Compute x := A^(-1) b, i.e., solve A x = b
void sbp_sat::x2::solve(
  KSP &A, std::vector<petsc_vector> &b, std::vector<petsc_vector> &x) {
    for (std::size_t i = 0; i < b.size(); ++i) {
        KSPSolve(A, b[i], x[i]);
    }
}

void sbp_sat::x2::fcompop(
  petsc_matrix &f, 
  petsc_matrix const &l, 
  petsc_matrix const &b,
  petsc_matrix const &h,
  real_t       const τ, 
  real_t       const β) {

  petsc_matrix t;
  MatMatMult(l, b, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &t);
  MatScale(t, β);
  MatAXPY(t, -τ, l, UNKNOWN_NONZERO_PATTERN);
  MatMatMult(t, h, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &f);
  finalize<fw>(f);
  destroy<fw>(t);
}

void sbp_sat::x2::make_f_subs(
  components const &sbp, 
  std::vector<std::vector<petsc_vector>> &f) {

  petsc_matrix f_n, f_s, f_e, f_w;

  // (-τ * LN + β * LN* BS_y) * H_x 
  fcompop(f_n, sbp.ln, sbp.bsy, sbp.hx, sbp.τ, sbp.β);
  // (-τ * LS + β * LS* BS_x) * H_x 
  fcompop(f_s, sbp.ls, sbp.bsy, sbp.hx, sbp.τ, sbp.β);
  // (-τ * LW + β * LE* BS_y) * H_y 
  fcompop(f_e, sbp.le, sbp.bsx, sbp.hy, sbp.τ, sbp.β);
  // (-τ * LW + β * LW* BS_y) * H_y 
  fcompop(f_w, sbp.lw, sbp.bsx, sbp.hy, sbp.τ, sbp.β);

  int ncols;
  int const *cols;
  const double *vals;

  f[0] = std::vector<petsc_vector>(sbp.n);
  f[1] = std::vector<petsc_vector>(sbp.n);
  f[2] = std::vector<petsc_vector>(sbp.n);
  f[3] = std::vector<petsc_vector>(sbp.n);

  for (std::size_t i = 0; i != sbp.n; ++i) {
    VecCreateSeq(PETSC_COMM_SELF, sbp.n * sbp.n, &f[0][i]);
    MatGetRow(f_e, i, &ncols, &cols, &vals);
    VecSetValues(f[0][i], ncols, cols, vals, ADD_VALUES);
    finalize<fw>(f[0][i]);

    VecCreateSeq(PETSC_COMM_SELF, sbp.n * sbp.n, &f[1][i]);
    MatGetRow(f_w, i, &ncols, &cols, &vals);
    VecSetValues(f[1][i], ncols, cols, vals, ADD_VALUES);
    finalize<fw>(f[1][i]);
  
    VecCreateSeq(PETSC_COMM_SELF, sbp.n * sbp.n, &f[2][i]);
    MatGetRow(f_s, i, &ncols, &cols, &vals);
    VecSetValues(f[2][i], ncols, cols, vals, ADD_VALUES);
    finalize<fw>(f[2][i]);
  
    VecCreateSeq(PETSC_COMM_SELF, sbp.n * sbp.n, &f[3][i]);
    MatGetRow(f_n, i, &ncols, &cols, &vals);
    VecSetValues(f[3][i], ncols, cols, vals, ADD_VALUES);
    finalize<fw>(f[3][i]);
  }

  destroy<fw>(f_n);
  destroy<fw>(f_s);
  destroy<fw>(f_e);
  destroy<fw>(f_w);
}

void sbp_sat::x2::msolvef(
  std::vector<std::vector<petsc_vector>> &x,
  KSP *m,
  std::size_t size,
  std::vector<std::vector<petsc_vector>> &f) {

  std::size_t index = 0;
  for (std::size_t i = 0; i != size; ++i) { //  ------------ n blocks
    for (std::size_t j = 0; j != f.size(); ++j) { // ------- 4
      for (std::size_t k = 0; k != f[j].size(); ++k) { // -- n blocks
        std::cout << "solve with block " << i << " on "
          << " f slice " << j << ", " << k << " into x index " 
          << index << ", " << k << std::endl;
        KSPSolve(m[i], f[j][k], x[index][k]);
      }
      index += 1;
    }
  }
}
 




