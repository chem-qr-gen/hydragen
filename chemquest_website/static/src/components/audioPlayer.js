const AudioPlayer = {
  audio: new Audio(),
  isPlaying: false,
  load: function (url){
    this.audio = new Audio(url);
  },
  play: function () {
    this.audio.play();
    this.isPlaying = true;
  },
  pause: function () {
    this.audio.pause();
    this.isPlaying = false;
  },
  toggle: function () {
    if (this.audio.paused) {
      this.play();
    } else {
      this.pause();
    }
  },
  view: () => {}
};

export default AudioPlayer;

