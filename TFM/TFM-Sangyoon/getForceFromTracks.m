function [forceMaturing, slopeMaturing] = getForceFromTracks(tracksNA)
    forceMaturing = cell(1,numel(tracksNA));
    slopeMaturing = zeros(1,numel(tracksNA));
    for k=1:numel(tracksNA)
        if ~islogical(tracksNA(k).presence)
            tracksNA(k).presence=logical(tracksNA(k).presence);
        end
        forceMaturing{k} = tracksNA(k).forceMag(tracksNA(k).presence);
        t = 5/60*(tracksNA(k).iFrame(tracksNA(k).presence));
        f = tracksNA(k).forceMag(tracksNA(k).presence);
        [~,slopeMaturing(k)] = regression(t,f);
    end
end
