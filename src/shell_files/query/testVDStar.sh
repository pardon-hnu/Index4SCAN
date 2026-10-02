

nohup ../VDStar/build/DynamicStrClu_ours -graph /home/hnu/Disk0/ParDon/graph/youtube.txt > ./res/temp_query/youtube/vd.txt &
nohup ./main youtube GS-Index exp-query > ./res/temp_query/youtube/gs-index.txt &
nohup ./main youtube BOTBIN exp-query > ./res/temp_query/youtube/botbin.txt & 
nohup ./main youtube FOREST exp-query > ./res/temp_query/youtube/forest.txt &
nohup ./main youtube PPT exp-query > ./res/temp_query/youtube/ppt.txt &

nohup ../VDStar/build/DynamicStrClu_ours -graph /home/hnu/Disk0/ParDon/graph/skitter.txt > ./res/temp_query/skitter/vd.txt &
nohup ./main skitter GS-Index exp-query > ./res/temp_query/skitter/skitter-gs-index.txt &
nohup ./main skitter BOTBIN exp-query > ./res/temp_query/skitter/skitter-botbin.txt &
nohup ./main skitter FOREST exp-query > ./res/temp_query/skitter/skitter-forest.txt &
nohup ./main skitter PPT exp-query > ./res/temp_query/skitter/skitter-ppt.txt &


nohup ../VDStar/build/DynamicStrClu_ours -graph /home/hnu/Disk0/ParDon/graph/topcats.txt > ./res/temp_query/topcats/vd.txt &
nohup ./main topcats GS-Index exp-query > ./res/temp_query/topcats/topcats-gs-index.txt &
nohup ./main topcats BOTBIN exp-query > ./res/temp_query/topcats/topcats-botbin.txt &
nohup ./main topcats FOREST exp-query > ./res/temp_query/topcats/topcats-forest.txt &
nohup ./main topcats PPT exp-query > ./res/temp_query/topcats/topcats-ppt.txt &

nohup ../VDStar/build/DynamicStrClu_ours -graph /home/hnu/Disk0/ParDon/graph/pokec.txt > ./res/temp_query/pokec/pokec-vd.txt &
nohup ./main pokec GS-Index exp-query > ./res/temp_query/pokec/pokec-gs-index.txt &
nohup ./main pokec BOTBIN exp-query > ./res/temp_query/pokec/pokec-botbin.txt &
nohup ./main pokec FOREST exp-query > ./res/temp_query/pokec/pokec-forest.txt &
nohup ./main pokec PPT exp-query > ./res/temp_query/pokec/pokec-ppt.txt &