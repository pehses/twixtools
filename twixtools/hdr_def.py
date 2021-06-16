import numpy as np

align_to_nbytes = 512
max_raidfile_entries = 64

# 8B
MrParcRaidFileHeader = [("hdSize_", "<u4"),     # ID???
                        ("count_", "<u4")]      # # of meas

# 152B
MrParcRaidFileEntry = [("measId_", "<u4"),      # MeasID
                       ("fileId_", "<u4"),      # FileID
                       ("off_", "<u8"),         # Measurement offset
                       ("len_", "<u8"),         # Measurement length
                       ("patName_", "<S64"),    # Patient name
                       ("protName_", "<S64")]   # Protocol name

# 9736B
MultiRaidFileHeader = [("hdr", MrParcRaidFileHeader),
                       ("entry", (MrParcRaidFileEntry, 64))]
MultiRaidFileHeader = np.dtype(MultiRaidFileHeader)
MrParcRaidFileHeader = np.dtype(MrParcRaidFileHeader)
MrParcRaidFileEntry = np.dtype(MrParcRaidFileEntry)


SingleMeasInit = [("hdr_len", "<u4"),
                  ("unknown", "<u4")]  # usually value 4
SingleMeasInit = np.dtype(SingleMeasInit)
